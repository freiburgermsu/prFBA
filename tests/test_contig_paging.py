"""Offline tests of the paged BV-BRC contig fetch (re-score decision D6, 2026-11-15).

The MAG-arm contig fetchers used to send ONE genome_sequence request with limit(50000); the server
caps a response at 25,000 rows, so 1733.9197 / 1733.9535 (random panel) were written truncated.
prefetch_gene_families.fetch_genome_sequence_rows now pages limit(25000,start) sorted by
sequence_id and raises unless the rows equal the Content-Range total. No test touches the network:
prefetch_gene_families._get is replaced by a fake BV-BRC serving pages from an in-memory list.
"""
import importlib.util
import os
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import prefetch_gene_families as PF  # noqa: E402

LIMIT_RE = re.compile(r"&limit\((\d+),(\d+)\)")
FETCHERS = ["mag_validation/fetch_contigs.py", "mag_validation_km/fetch_contigs_km.py",
            "mag_validation_rnd/fetch_contigs_rnd.py", "mag_validation_rnd/fetch_contigs_rnd_full.py"]


def contigs(n, gid="1733.1"):
    return [{"sequence_id": f"{gid}.con.{i:06d}", "accession": f"{gid}.con.{i:06d}", "sequence": "ACGT"}
            for i in range(n)]


class FakeBVBRC:
    """Serves limit(n,start) pages of `rows`; cap = the server's max rows per response;
    ignore_start = a server that ignores the offset (returns the first page every time)."""

    def __init__(self, rows, cap=25000, total=None, ignore_start=False, no_range=False):
        self.rows, self.cap, self.ignore_start, self.no_range = rows, cap, ignore_start, no_range
        self.total = len(rows) if total is None else total
        self.urls = []

    def __call__(self, url, timeout=60, retries=3):
        self.urls.append(url)
        assert "genome_sequence" in url and "sort(+sequence_id)" in url and "select(sequence_id," in url
        n, start = map(int, LIMIT_RE.search(url).groups())
        if self.ignore_start:
            start = 0
        page = self.rows[start:start + min(n, self.cap)]
        return page, (None if self.no_range else self.total)


@pytest.mark.parametrize("n", [0, 1, 24999, 25000, 25001, 36062, 75000])
def test_pages_until_complete(monkeypatch, n):
    fake = FakeBVBRC(contigs(n))
    monkeypatch.setattr(PF, "_get", fake)
    rows = PF.fetch_genome_sequence_rows("1733.1")
    assert rows == contigs(n)
    assert len(fake.urls) == max(1, -(-n // PF.PAGE_SIZE))      # ceil(n / 25,000) requests


def test_the_old_cap_is_detected(monkeypatch):
    # a server that caps at 25,000 and ignores the offset would repeat rows: refuse
    fake = FakeBVBRC(contigs(31169), ignore_start=True)
    monkeypatch.setattr(PF, "_get", fake)
    with pytest.raises(PF.IncompleteFetch):
        PF.fetch_genome_sequence_rows("1733.9197")


def test_short_page_before_total_raises(monkeypatch):
    fake = FakeBVBRC(contigs(30000), total=31169)
    monkeypatch.setattr(PF, "_get", fake)
    with pytest.raises(PF.IncompleteFetch, match="server reports 31169"):
        PF.fetch_genome_sequence_rows("1733.9197")


def test_missing_content_range_raises(monkeypatch):
    monkeypatch.setattr(PF, "_get", FakeBVBRC(contigs(10), no_range=True))
    with pytest.raises(PF.IncompleteFetch, match="no Content-Range"):
        PF.fetch_genome_sequence_rows("1733.1")


def test_total_changing_between_pages_raises(monkeypatch):
    fake = FakeBVBRC(contigs(30000))
    calls = {"n": 0}

    def flaky(url, timeout=60, retries=3):
        page, total = fake(url, timeout, retries)
        calls["n"] += 1
        return page, total + (calls["n"] > 1)
    monkeypatch.setattr(PF, "_get", flaky)
    with pytest.raises(PF.IncompleteFetch, match="changed between pages"):
        PF.fetch_genome_sequence_rows("1733.1")


@pytest.mark.parametrize("path", FETCHERS)
def test_every_contig_fetcher_pages(monkeypatch, tmp_path, path):
    spec = importlib.util.spec_from_file_location("fetcher_" + os.path.basename(path)[:-3],
                                                  os.path.join(ROOT, path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert "&limit(50000)" not in open(os.path.join(ROOT, path)).read()   # the old unpaged URL
    fake = FakeBVBRC(contigs(36001, "1733.9535"))
    monkeypatch.setattr(mod.PF, "_get", fake)
    rows = mod.fetch("1733.9535")
    assert len(rows) == 36001 and len(fake.urls) == 2
    monkeypatch.setattr(mod, "OUT", str(tmp_path))
    mod.write("1733.9535", rows)
    fna = (tmp_path / "1733.9535.fna").read_text()
    assert fna.count(">") == 36001 and fna.startswith(">1733.9535|1733.9535.con.000000\nACGT\n")
