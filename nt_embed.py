"""
nt_embed.py — shared Nucleotide-Transformer-v2 embedder for 16S rRNA sequences.

Model:   InstaDeepAI/nucleotide-transformer-v2-500m-multi-species
         498M params, 1024-d hidden, 6-mer tokenizer (~6 bp / token, 2048-token ctx).
Embedding:
         mean-pool of the final hidden layer over real (non-pad) tokens, then L2-normalize.
         L2-normalized vectors -> cosine similarity == dot product.

Why this recipe (Blackwell / RTX 5070 Ti, torch cu13, transformers 4.44.2):
  * NT-v2 uses a *gated* FFN and ships custom remote code; it only loads correctly
    through `trust_remote_code=True` on transformers 4.x (5.x forces the in-library
    ESM class, whose plain FFN shape-mismatches the gated checkpoint).
  * The remote attention upcasts the softmax to fp32 while values stay in the storage
    dtype, so loading weights directly in bf16 raises a matmul dtype clash. We therefore
    keep *fp32 master weights* and wrap the forward pass in `torch.autocast(bf16)`, which
    reconciles operand dtypes and still runs the matmuls on bf16 tensor cores.

This module is imported by embed_16s.py (bulk embedding) and query.py (reference lookup)
so both use byte-for-byte the same encoder.
"""
from __future__ import annotations
import torch
import torch.nn.functional as F

MODEL_ID = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
EMBED_DIM = 1024


def load_model(model_id: str = MODEL_ID, device: str = "cuda"):
    """Load tokenizer + model (fp32 weights, eval mode) ready for autocast inference."""
    from transformers import AutoTokenizer, AutoModelForMaskedLM
    tok = AutoTokenizer.from_pretrained(model_id)
    model = AutoModelForMaskedLM.from_pretrained(model_id, trust_remote_code=True)
    model = model.to(device).eval()
    return tok, model


@torch.no_grad()
def embed_batch(seqs, tok, model, device: str = "cuda", max_tokens: int = 512) -> torch.Tensor:
    """
    Embed a list of DNA strings -> (len(seqs), 1024) L2-normalized float32 tensor on `device`.

    max_tokens caps the tokenized length (512 tokens ~ 3 kb): full-length 16S (~1.5 kb,
    ~260 tokens) is never truncated; only atypical multi-kb mis-exports are clipped.
    """
    enc = tok.batch_encode_plus(
        list(seqs), return_tensors="pt", padding=True,
        truncation=True, max_length=max_tokens,
    )
    enc = {k: v.to(device) for k, v in enc.items()}
    with torch.autocast(device_type=device.split(":")[0], dtype=torch.bfloat16):
        out = model(**enc, output_hidden_states=True)
    h = out.hidden_states[-1].float()                      # (B, T, 1024)
    mask = enc["attention_mask"].unsqueeze(-1).float()     # (B, T, 1)
    summed = (h * mask).sum(dim=1)
    counts = mask.sum(dim=1).clamp(min=1.0)
    emb = summed / counts                                  # masked mean pool
    return F.normalize(emb, p=2, dim=1)                    # unit vectors -> cosine = dot
