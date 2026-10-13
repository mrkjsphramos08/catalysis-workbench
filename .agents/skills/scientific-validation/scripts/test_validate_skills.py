"""Validator regression tests using isolated temporary directories only."""

from pathlib import Path
import tempfile
import unittest

from validate_skills import validate_links, validate_repository, validate_skill


class SkillValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.skill = self.root / "example-skill"
        self.skill.mkdir()

    def write_manifest(self, frontmatter: str, body: str = "Read the relevant source before editing.") -> None:
        (self.skill / "SKILL.md").write_text(f"---\n{frontmatter}\n---\n{body}\n", encoding="utf-8")

    def test_valid_yaml_and_directory(self) -> None:
        self.write_manifest("name: example-skill\ndescription: |\n  Reviews experimental measurements when auditing data.")
        self.assertEqual(validate_skill(self.skill), [])

    def test_rejects_invalid_frontmatter(self) -> None:
        cases = [
            "name: wrong-name\ndescription: Valid description.",
            "name: example-skill\ndescription: ''",
            "name: example-skill\ndescription: [not, a, string]",
            "name: example-skill\ndescription: valid\nname: example-skill",
            "name: example-skill\ndescription: valid\nvendor-trigger: true",
            "name: example-skill\ndescription: [unterminated",
            "- list\n- instead-of-mapping",
        ]
        for case in cases:
            with self.subTest(frontmatter=case):
                self.write_manifest(case)
                self.assertTrue(validate_skill(self.skill))

    def test_missing_manifest_and_empty_body(self) -> None:
        self.assertTrue(validate_skill(self.skill))
        self.write_manifest("name: example-skill\ndescription: Valid description.", body="")
        self.assertTrue(validate_skill(self.skill))

    def test_local_links_and_encoded_filename(self) -> None:
        reference = self.skill / "reference file.md"
        reference.write_text("Reference", encoding="utf-8")
        document = self.skill / "SKILL.md"
        document.write_text("[valid](reference%20file.md) [external](https://example.org) [heading](#local)", encoding="utf-8")
        self.assertEqual(validate_links(document, self.root), [])
        document.write_text("[missing](missing.md) [escape](../../outside.md)", encoding="utf-8")
        self.assertEqual(len(validate_links(document, self.root)), 2)

    def test_required_collection_is_not_vacuously_valid(self) -> None:
        errors, count, _ = validate_repository(self.root)
        self.assertEqual(count, 0)
        self.assertEqual(sum(error.startswith("Missing required skill:") for error in errors), 7)


if __name__ == "__main__":
    unittest.main()
