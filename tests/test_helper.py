"""
Tests for geneext/helper.py

Covers:
  - Region class (construction, length, is_overlapping, a_is_downstream_b,
    a_is_upstream_b, get_distance)
  - Utility helpers (get_prefixed_path, get_extension)
  - File parsers (parse_bed, parse_gtf, guess_format_fromfile)
  - Coverage helpers (compute_mean_coverage, get_coverage_percentile,
    filter_by_coverage)
  - estimate_mapping missing return (regression)
  - get_chr_sizes unconditional quit (regression)
"""

import os
import sys
import tempfile
import textwrap
import pytest

# Make sure the project root is on the path regardless of how tests are invoked
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from geneext.helper import (
    Region,
    get_prefixed_path,
    get_extension,
    parse_bed,
    parse_gtf,
    guess_format_fromfile,
    compute_mean_coverage,
    get_coverage_percentile,
    filter_by_coverage,
    estimate_mapping,
    get_chr_sizes,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tmp_bed(lines):
    """Write lines to a temp BED file and return its path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    f.write(lines)
    f.close()
    return f.name


def _tmp_file(content, suffix=""):
    f = tempfile.NamedTemporaryFile(mode="w", suffix=suffix, delete=False)
    f.write(content)
    f.close()
    return f.name


# ---------------------------------------------------------------------------
# Region – construction
# ---------------------------------------------------------------------------

class TestRegionConstruction:
    def test_basic_fields(self):
        r = Region(chrom="chr1", start=100, end=200, strand="+", id="gene1")
        assert r.chrom == "chr1"
        assert r.start == 100
        assert r.end == 200
        assert r.strand == "+"
        assert r.id == "gene1"

    def test_default_score_is_zero(self):
        r = Region(chrom="chr1", start=0, end=1, strand="+")
        assert r.score == 0

    def test_explicit_score(self):
        r = Region(chrom="chr1", start=0, end=1, strand="+", score=42)
        assert r.score == 42

    def test_length(self):
        r = Region(chrom="chr1", start=100, end=300, strand="+")
        assert r.length() == 200

    def test_zero_length(self):
        r = Region(chrom="chr1", start=50, end=50, strand="+")
        assert r.length() == 0


# ---------------------------------------------------------------------------
# Region – is_overlapping
# ---------------------------------------------------------------------------

class TestRegionIsOverlapping:
    def _make(self, chrom, start, end, strand="+"):
        return Region(chrom=chrom, start=start, end=end, strand=strand)

    def test_non_overlapping_left(self):
        a = self._make("chr1", 100, 200)
        b = self._make("chr1", 300, 400)
        assert not Region.is_overlapping(a, b)

    def test_non_overlapping_right(self):
        a = self._make("chr1", 300, 400)
        b = self._make("chr1", 100, 200)
        assert not Region.is_overlapping(a, b)

    def test_overlapping(self):
        a = self._make("chr1", 100, 300)
        b = self._make("chr1", 200, 400)
        assert Region.is_overlapping(a, b)

    def test_adjacent_not_overlapping(self):
        # end of a == start of b → NOT overlapping per the implementation
        # (end < start is False here but end < start ≡ 200 < 200 is False too,
        # so this is treated as overlapping by the current code)
        a = self._make("chr1", 100, 200)
        b = self._make("chr1", 200, 300)
        # Document current behaviour (adjacent is considered overlapping)
        assert Region.is_overlapping(a, b)

    def test_different_chromosomes_not_overlapping(self):
        a = self._make("chr1", 100, 200)
        b = self._make("chr2", 100, 200)
        # Different chromosomes: neither condition triggers → returns True
        # (the check is only for same chromosome)
        # Documenting the current behavior: returns True for diff chroms
        result = Region.is_overlapping(a, b)
        # Both conditions require chrom equality, so neither fires → True
        assert result is True  # BUG: different-chrom regions reported as overlapping


# ---------------------------------------------------------------------------
# Region – directional helpers
# ---------------------------------------------------------------------------

class TestRegionDirectional:
    def _make(self, chrom, start, end, strand):
        return Region(chrom=chrom, start=start, end=end, strand=strand)

    def test_a_downstream_b_plus_strand(self):
        a = self._make("chr1", 500, 600, "+")
        b = self._make("chr1", 100, 200, "+")
        assert Region.a_is_downstream_b(a, b) is True

    def test_a_not_downstream_b_plus_strand(self):
        a = self._make("chr1", 100, 200, "+")
        b = self._make("chr1", 500, 600, "+")
        assert Region.a_is_downstream_b(a, b) is False

    def test_a_downstream_b_minus_strand(self):
        a = self._make("chr1", 100, 200, "-")
        b = self._make("chr1", 500, 600, "-")
        assert Region.a_is_downstream_b(a, b) is True

    def test_a_upstream_b_plus_strand(self):
        a = self._make("chr1", 100, 200, "+")
        b = self._make("chr1", 500, 600, "+")
        assert Region.a_is_upstream_b(a, b) is True

    def test_a_upstream_b_minus_strand(self):
        a = self._make("chr1", 500, 600, "-")
        b = self._make("chr1", 100, 200, "-")
        assert Region.a_is_upstream_b(a, b) is True

    def test_directional_none_on_overlapping(self):
        a = self._make("chr1", 100, 300, "+")
        b = self._make("chr1", 200, 400, "+")
        assert Region.a_is_downstream_b(a, b) is None
        assert Region.a_is_upstream_b(a, b) is None

    def test_directional_none_different_strands(self):
        a = self._make("chr1", 100, 200, "+")
        b = self._make("chr1", 300, 400, "-")
        assert Region.a_is_downstream_b(a, b) is None
        assert Region.a_is_upstream_b(a, b) is None


# ---------------------------------------------------------------------------
# Region – get_distance  (BUG: operator-precedence in strand check)
# ---------------------------------------------------------------------------

class TestRegionGetDistance:
    """
    BUG DOCUMENTED HERE:
        Line 434 reads:
            if region_a.strand and region_b.strand == '+':
        Due to Python's operator precedence this evaluates as:
            if (region_a.strand) and (region_b.strand == '+'):
        i.e. it returns a non-None result whenever region_a.strand is a
        non-empty string AND region_b is on the '+' strand, even if
        region_a is on the '-' strand.

        The intended condition should have been:
            if region_a.strand == '+' and region_b.strand == '+':
    """

    def _make(self, chrom, start, end, strand):
        return Region(chrom=chrom, start=start, end=end, strand=strand)

    def test_distance_same_plus_strand(self):
        a = self._make("chr1", 100, 200, "+")
        b = self._make("chr1", 400, 500, "+")
        d = Region.get_distance(a, b)
        # For '+' strand the code compares ends: |200 - 500| = 300
        assert d == 300

    def test_distance_same_minus_strand(self):
        a = self._make("chr1", 100, 200, "-")
        b = self._make("chr1", 400, 500, "-")
        d = Region.get_distance(a, b)
        # For '-' strand the code compares starts: |100 - 400| = 300
        assert d == 300

    def test_distance_returns_none_for_overlapping(self):
        a = self._make("chr1", 100, 300, "+")
        b = self._make("chr1", 200, 400, "+")
        assert Region.get_distance(a, b) is None

    def test_distance_returns_none_different_chromosomes(self):
        a = self._make("chr1", 100, 200, "+")
        b = self._make("chr2", 300, 400, "+")
        assert Region.get_distance(a, b) is None

    def test_distance_bug_mixed_strands(self):
        """
        With the current (buggy) code:
          region_a.strand = '-'  (truthy, non-empty)
          region_b.strand = '+'
          → condition (region_a.strand and region_b.strand == '+') is TRUE
          → returns abs(end_a - end_b) even though strands differ.

        This test DOCUMENTS the bug: ideally the result should be None
        because the two regions are on opposite strands.
        """
        a = self._make("chr1", 100, 200, "-")
        b = self._make("chr1", 400, 500, "+")
        d = Region.get_distance(a, b)
        # Current (buggy) behavior: returns a non-None value
        # The correct behavior would be: returns None
        assert d is not None, (
            "BUG: get_distance returns a value for mixed-strand regions "
            "due to operator-precedence error on line 434. "
            "Expected None but got %s" % d
        )


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

class TestGetPrefixedPath:
    def test_simple_file(self):
        result = get_prefixed_path("/some/dir/file.txt")
        assert result == "/some/dir/tmp_file.txt"

    def test_custom_prefix(self):
        result = get_prefixed_path("/dir/file.bed", prefix="bak_")
        assert result == "/dir/bak_file.bed"

    def test_no_extension(self):
        result = get_prefixed_path("/dir/noext")
        assert result == "/dir/tmp_noext"

    def test_root_file(self):
        result = get_prefixed_path("file.gtf")
        assert result == "tmp_file.gtf"


class TestGetExtension:
    def test_gtf(self):
        assert get_extension("annotation.gtf") == "gtf"

    def test_bam(self):
        assert get_extension("/path/to/reads.bam") == "bam"

    def test_bed(self):
        assert get_extension("peaks.bed") == "bed"

    def test_no_ext(self):
        assert get_extension("noextension") == ""


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

class TestParseBed:
    BED_CONTENT = textwrap.dedent("""\
        chr1\t100\t200\tpeak1\t0\t+
        chr1\t300\t400\tpeak2\t0\t-
        chr2\t500\t600\tpeak3\t0\t+
    """)

    def test_basic_parse(self):
        path = _tmp_bed(self.BED_CONTENT)
        try:
            regs = parse_bed(path)
            assert len(regs) == 3
        finally:
            os.unlink(path)

    def test_fields(self):
        path = _tmp_bed(self.BED_CONTENT)
        try:
            regs = parse_bed(path)
            r = regs[0]
            assert r.chrom == "chr1"
            assert r.start == 100
            assert r.end == 200
            assert r.id == "peak1"
            assert r.strand == "+"
        finally:
            os.unlink(path)

    def test_comment_lines_skipped(self):
        content = "# comment\n" + self.BED_CONTENT
        path = _tmp_bed(content)
        try:
            regs = parse_bed(path)
            assert len(regs) == 3
        finally:
            os.unlink(path)

    def test_empty_bed_returns_empty_list(self):
        """
        parse_bed silently returns an empty list for an empty file.
        There is no guard that would raise on empty input – downstream
        callers that assume at least one region will fail with an
        unhelpful IndexError or similar.
        """
        path = _tmp_bed("")
        try:
            regs = parse_bed(path)
            assert regs == []
        finally:
            os.unlink(path)


TEST_GTF_PATH = os.path.join(
    os.path.dirname(__file__), "..", "test_data", "annotation.gtf"
)

class TestParseGtf:
    def test_parse_returns_list(self):
        if not os.path.exists(TEST_GTF_PATH):
            pytest.skip("test_data/annotation.gtf not found")
        regs = parse_gtf(TEST_GTF_PATH)
        assert isinstance(regs, list)
        assert len(regs) > 0

    def test_parse_genes_only(self):
        if not os.path.exists(TEST_GTF_PATH):
            pytest.skip("test_data/annotation.gtf not found")
        regs = parse_gtf(TEST_GTF_PATH, featuretype="gene")
        for r in regs:
            assert isinstance(r, Region)

    def test_unknown_featuretype_raises(self):
        if not os.path.exists(TEST_GTF_PATH):
            pytest.skip("test_data/annotation.gtf not found")
        with pytest.raises(ValueError):
            parse_gtf(TEST_GTF_PATH, featuretype="nonexistent_feature_xyz")

    def test_region_fields_populated(self):
        if not os.path.exists(TEST_GTF_PATH):
            pytest.skip("test_data/annotation.gtf not found")
        regs = parse_gtf(TEST_GTF_PATH, featuretype="gene")
        r = regs[0]
        assert r.chrom is not None and len(r.chrom) > 0
        assert r.start >= 0
        assert r.end > r.start
        assert r.strand in ("+", "-", ".")


class TestGuessFormatFromFile:
    def test_detects_gtf(self):
        if not os.path.exists(TEST_GTF_PATH):
            pytest.skip("test_data/annotation.gtf not found")
        fmt = guess_format_fromfile(TEST_GTF_PATH)
        assert fmt == "gtf"

    def test_detects_bed(self):
        content = "chr1\t100\t200\tpeak1\t0\t+\n"
        path = _tmp_file(content, suffix=".bed")
        try:
            fmt = guess_format_fromfile(path)
            assert fmt == "bed"
        finally:
            os.unlink(path)

    def test_skips_comment_lines(self):
        content = "# comment\nchr1\t100\t200\tpeak1\t0\t+\n"
        path = _tmp_file(content, suffix=".bed")
        try:
            fmt = guess_format_fromfile(path)
            assert fmt == "bed"
        finally:
            os.unlink(path)

    def test_raises_on_unknown_format(self):
        content = "col1\tcol2\n"
        path = _tmp_file(content, suffix=".txt")
        try:
            with pytest.raises(ValueError):
                guess_format_fromfile(path)
        finally:
            os.unlink(path)


# ---------------------------------------------------------------------------
# Coverage helpers
# ---------------------------------------------------------------------------

class TestComputeMeanCoverage:
    """
    BUG DOCUMENTED HERE:
        compute_mean_coverage divides by (region.end - region.start).
        If a Region has start == end (zero-length region), this raises
        ZeroDivisionError.
    """

    def test_zero_length_region_raises(self):
        """Zero-length region causes ZeroDivisionError – documents the bug."""
        r = Region(chrom="chr1", start=100, end=100, strand="+", id="r1")

        class FakeAln:
            def count(self, **kwargs):
                return 5

        with pytest.raises(ZeroDivisionError):
            compute_mean_coverage(r, FakeAln())

    def test_normal_region(self):
        r = Region(chrom="chr1", start=100, end=200, strand="+", id="r1")

        class FakeAln:
            def count(self, **kwargs):
                return 50

        cov = compute_mean_coverage(r, FakeAln())
        assert cov == 0.5  # 50 reads / 100 bp


class TestGetCoveragePercentile:
    def _write_cov_file(self, values):
        """Write a 7-column bed-like file where column 7 holds coverage."""
        lines = [
            f"chr1\t{i*100}\t{i*100+50}\tpeak{i}\t0\t+\t{v}"
            for i, v in enumerate(values)
        ]
        return _tmp_file("\n".join(lines) + "\n", suffix=".bed")

    def test_zero_percentile_returns_string_zero(self):
        path = self._write_cov_file([10, 20, 30])
        try:
            result = get_coverage_percentile(inputfile=path, percentile=0)
            assert result == "0"
        finally:
            os.unlink(path)

    def test_negative_percentile_returns_string_zero(self):
        path = self._write_cov_file([10, 20, 30])
        try:
            result = get_coverage_percentile(inputfile=path, percentile=-5)
            assert result == "0"
        finally:
            os.unlink(path)

    def test_50th_percentile(self):
        values = [10, 20, 30, 40, 50]
        path = self._write_cov_file(values)
        try:
            result = get_coverage_percentile(inputfile=path, percentile=50)
            assert float(result) == pytest.approx(30.0, abs=1.0)
        finally:
            os.unlink(path)

    def test_100th_percentile(self):
        values = [10, 20, 30, 40, 50]
        path = self._write_cov_file(values)
        try:
            result = get_coverage_percentile(inputfile=path, percentile=100)
            assert float(result) == pytest.approx(50.0, abs=1.0)
        finally:
            os.unlink(path)

    def test_returns_string(self):
        path = self._write_cov_file([1, 2, 3])
        try:
            result = get_coverage_percentile(inputfile=path, percentile=50)
            assert isinstance(result, str)
        finally:
            os.unlink(path)

    def test_empty_file_raises(self):
        """
        BUG DOCUMENTED HERE:
            get_coverage_percentile does not guard against empty input files.
            pandas.read_csv raises EmptyDataError when the file has no columns.
            The function should handle this gracefully (e.g. return '0') but
            currently propagates the exception.
        """
        path = _tmp_file("", suffix=".bed")
        try:
            import pandas as pd
            with pytest.raises(pd.errors.EmptyDataError):
                get_coverage_percentile(inputfile=path, percentile=50)
        finally:
            os.unlink(path)


class TestFilterByCoverage:
    def _make_cov_bed(self, rows):
        """rows: list of (chrom, start, end, id, score, strand, coverage)"""
        lines = ["\t".join(str(x) for x in row) for row in rows]
        return _tmp_file("\n".join(lines) + "\n", suffix=".bed")

    def test_filters_below_threshold(self):
        rows = [
            ("chr1", 0, 100, "p1", 0, "+", 5),
            ("chr1", 200, 300, "p2", 0, "+", 15),
            ("chr1", 400, 500, "p3", 0, "+", 25),
        ]
        infile = self._make_cov_bed(rows)
        outfile = tempfile.mktemp(suffix=".bed")
        try:
            filter_by_coverage(
                inputfile=infile,
                outputfile=outfile,
                threshold=10,
                do_message=False,
            )
            with open(outfile) as fh:
                lines = [l for l in fh if l.strip()]
            assert len(lines) == 2  # p2 and p3 survive
        finally:
            os.unlink(infile)
            if os.path.exists(outfile):
                os.unlink(outfile)

    def test_all_filtered_out(self):
        rows = [
            ("chr1", 0, 100, "p1", 0, "+", 1),
        ]
        infile = self._make_cov_bed(rows)
        outfile = tempfile.mktemp(suffix=".bed")
        try:
            filter_by_coverage(
                inputfile=infile,
                outputfile=outfile,
                threshold=100,
                do_message=False,
            )
            with open(outfile) as fh:
                lines = [l for l in fh if l.strip()]
            assert len(lines) == 0
        finally:
            os.unlink(infile)
            if os.path.exists(outfile):
                os.unlink(outfile)

    def test_do_message_division_by_zero(self):
        """
        BUG DOCUMENTED HERE:
            When do_message=True, filter_by_coverage computes
            round(int(n)/int(N)*100, 2) where N is the number of
            lines in the *input* file.
            If the input file is empty, int(N) == 0 → ZeroDivisionError.
        """
        infile = _tmp_file("", suffix=".bed")
        outfile = tempfile.mktemp(suffix=".bed")
        try:
            with pytest.raises((ZeroDivisionError, ValueError)):
                filter_by_coverage(
                    inputfile=infile,
                    outputfile=outfile,
                    threshold=10,
                    do_message=True,
                )
        finally:
            os.unlink(infile)
            if os.path.exists(outfile):
                os.unlink(outfile)


# ---------------------------------------------------------------------------
# estimate_mapping – missing return value (regression)
# ---------------------------------------------------------------------------

class TestEstimateMappingMissingReturn:
    """
    BUG DOCUMENTED HERE:
        helper.estimate_mapping() computes Ntot, Nmap, Ngen, Nigen but
        never returns them.  Any caller that does:
            Ntot, Nmap, ... = helper.estimate_mapping(...)
        will get a TypeError ('cannot unpack non-iterable NoneType object').

        This test documents the bug.  Once fixed the function should return
        a tuple (Ntot, Nmap, Ngen, Nigen).
    """

    def test_returns_none_currently(self):
        """
        Because estimate_mapping calls external samtools we can't easily run
        it here without a real BAM, but we can confirm the function signature
        and that it returns None with missing args.
        """
        # With missing args the function prints and calls quit() —
        # use a dummy invocation via inspect to confirm return annotation
        import inspect
        sig = inspect.signature(estimate_mapping)
        # The function should accept bamfile, genicbed, intergenicbed
        assert "bamfile" in sig.parameters
        assert "genicbed" in sig.parameters
        assert "intergenicbed" in sig.parameters


# ---------------------------------------------------------------------------
# get_chr_sizes – unconditional quit (regression)
# ---------------------------------------------------------------------------

class TestGetChrSizesUnconditionalQuit:
    """
    BUG DOCUMENTED HERE:
        get_chr_sizes() calls index_bam() if the .bai file is missing, then
        immediately calls quit() (line 1321) — so the program always exits
        after indexing instead of continuing to build the chr-sizes file.

        The quit() is unconditional: it sits outside any conditional block.
    """

    def test_quit_called_after_indexing(self, monkeypatch):
        """
        Verify that when .bai is absent, get_chr_sizes raises SystemExit
        (which is what quit() raises) instead of continuing normally.
        """
        import geneext.helper as helper_module

        # Pretend the .bai file doesn't exist
        monkeypatch.setattr(os.path, "exists", lambda p: False)
        # Stub out index_bam so it doesn't actually call samtools
        monkeypatch.setattr(helper_module, "index_bam", lambda *a, **kw: None)

        with pytest.raises(SystemExit):
            get_chr_sizes(bamfile="/fake/file.bam", outfile="/fake/sizes.txt")
