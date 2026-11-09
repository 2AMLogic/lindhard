"""Parser tests on a synthetic text fixture in the `pdftotext -layout` shape of
the srim.org citation files. Run: python3 -m unittest test_fetch_bibliography
(no network, no PDFs; every name and title below is invented)."""
import unittest

import fetch_bibliography as fb

FIXTURE = """\
Stopping for Ion =H                                              , Target = Xx
  Pub.                                                                                          Citation
  Year        Authors, Title, Journal Citation and Comments                                      Numb
           Alpha, A. B.                                                                        1962-Alph
1962 'A Synthetic Title for Testing'
     Phys. Rev. Letters, 10, 399-400 (1962)
                                                                                                 0513
Comment : S. 70 keV H -> Xx
           Beta, C.    Gamma, D. E.     Delta, F.                                              1970-Beta2
1970      'A Title That Wraps Onto a Second Line in the
          Layout of the File'
                                                                                                 0171


          Nucl. Inst. Methods, 80, 1-9 (1970)
Comment : S, dS. 2 MeV H -> Xx (Cryst.), a comment that
          continues on a second line
           Epsilon, G.                                                                         1975-Epsi
1975 'No Journal Line Follows This Title'
Comment : S. 1 MeV H -> Xx

(C) 2011 SRIM.org                                              Experimental Stopping and Range Papers
                                 Page 1 of 1
"""

COMPOUND = """\
    Citations for Target : Xx2O3
           Zeta, H.        Eta, I.                                                    1962-Zeta
1962       'A Compound Entry'                                                                  
           Can. J. Phys., 40, 1749-64 (1962)
 Comment : S. 4-30 keV H -> Xx2O3
"""


class ParseTests(unittest.TestCase):
    def parsed(self, text):
        return [fb.parse_block(b) for b in fb.split_entries(text)]

    def test_parsed_entry(self):
        a = self.parsed(FIXTURE)[0]
        self.assertEqual(a["citation_key"], "1962-Alph")
        self.assertEqual(a["year"], "1962")
        self.assertEqual(a["authors"], "Alpha, A. B.")
        self.assertEqual(a["title"], "A Synthetic Title for Testing")
        self.assertEqual(a["journal"], "Phys. Rev. Letters, 10, 399-400 (1962)")
        self.assertEqual(a["old_ref"], "0513")
        self.assertEqual(a["comment"], "S. 70 keV H -> Xx")

    def test_multi_line_entry(self):
        b = self.parsed(FIXTURE)[1]
        self.assertEqual(b["citation_key"], "1970-Beta2")
        self.assertEqual(b["authors"], "Beta, C. Gamma, D. E. Delta, F.")
        self.assertEqual(b["title"], "A Title That Wraps Onto a Second Line in the Layout of the File")
        self.assertEqual(b["journal"], "Nucl. Inst. Methods, 80, 1-9 (1970)")
        self.assertEqual(b["old_ref"], "0171")
        self.assertEqual(b["comment"], "S, dS. 2 MeV H -> Xx (Cryst.), a comment that continues on a second line")

    def test_unparsable_entry_is_kept_not_dropped(self):
        parsed = self.parsed(FIXTURE)
        self.assertEqual(len(parsed), 3)
        self.assertIsNone(parsed[2])  # no journal line
        self.assertEqual(fb.count_entries(FIXTURE), (3, 3))

    def test_compound_entry_without_old_ref(self):
        (c,) = self.parsed(COMPOUND)
        self.assertEqual(c["citation_key"], "1962-Zeta")
        self.assertEqual(c["old_ref"], "")
        rows = {}
        fb.merge(rows, c, target="Xx2O3")
        fb.merge(rows, c, pair="01-14")
        (r,) = rows.values()
        self.assertEqual((r["pairs"], r["target"]), (["01-14"], ["Xx2O3"]))

    def test_shared_entry_is_one_row_with_several_pairs(self):
        a = self.parsed(FIXTURE)[0]
        rows = {}
        for pair in ("01-14", "02-14", "01-14"):
            fb.merge(rows, dict(a), pair=pair)
        (r,) = rows.values()
        self.assertEqual(r["pairs"], ["01-14", "02-14"])


class DoiMatchTests(unittest.TestCase):
    def test_journal_abbreviation_matches_full_title(self):
        self.assertTrue(fb.journal_matches("Phys. Rev. Letters", ["Physical Review Letters"]))
        self.assertFalse(fb.journal_matches("Phys. Rev.", ["Physical Review Letters"]))
        self.assertFalse(fb.journal_matches("Phys. Letters", ["Physical Review Letters"]))

    def test_first_page_and_journal_name(self):
        ref = "Phys. Letters, 100A, 3, 137-140 (1984)"
        self.assertEqual(fb.first_page(ref), "137")
        self.assertEqual(fb.journal_name(ref), "Phys. Letters")

    def test_ambiguous_or_mismatched_candidates_leave_doi_empty(self):
        row = {"citation_key": "k", "title": "A Title", "year": "1962", "authors": "Alpha, A.",
               "journal": "Phys. Rev. Letters, 10, 399-400 (1962)"}
        good = {"DOI": "10.1/a", "title": ["A title"], "issued": {"date-parts": [[1962]]},
                "container-title": ["Physical Review Letters"], "page": "399-400"}
        for items, want in (
            ([good], "10.1/a"),
            ([good, dict(good, DOI="10.1/b")], ""),
            ([dict(good, page="401")], ""),
            ([dict(good, issued={"date-parts": [[1963]]})], ""),
            ([dict(good, title=["Another title"])], ""),
        ):
            cache = {"k|" + fb.norm(row["title"]): items}
            self.assertEqual(fb.crossref_doi(row, cache, [0]), want)


if __name__ == "__main__":
    unittest.main()
