"""Read-only checks for repository skills; not a scientific or activation validator.

Requires PyYAML in the selected interpreter. Checks inline Markdown file links,
not URL reachability, heading anchors, reference-style links, or API semantics.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

try:
    import yaml
except ImportError:
    raise SystemExit("PyYAML is required. Use the project interpreter or install PyYAML in your chosen environment.")


REQUIRED_SKILLS = {
    "scientific-python", "quantum-espresso", "electrochemistry",
    "materials-characterization", "scientific-validation",
    "publication-figures", "research-reproducibility",
}
STANDARD_FIELDS = {"name", "description", "license", "compatibility", "metadata", "allowed-tools"}
NAME = re.compile(r"[a-z]+(?:-[a-z]+)*\Z")  # This repository's letters-and-hyphens convention.
FRONTMATTER = re.compile(r"\A---\r?\n(.*?)\r?\n---(?:\r?\n|\Z)", re.DOTALL)
INLINE_LINK = re.compile(r"\[[^\]\n]*\]\((<[^>]+>|[^\s)]+)(?:\s+\"[^\"]*\")?\)")


class UniqueKeyLoader(yaml.SafeLoader):
    """Reject duplicate YAML fields rather than silently taking the last value."""


def _unique_mapping(loader: UniqueKeyLoader, node: yaml.MappingNode) -> dict:
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node)
        if not isinstance(key, str):
            raise ValueError("YAML mapping keys must be strings")
        if key in mapping:
            raise ValueError(f"Duplicate YAML key: {key}")
        mapping[key] = loader.construct_object(value_node)
    return mapping


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def validate_skill(folder: Path) -> list[str]:
    """Return structural problems for one skill, without executing its contents."""
    manifest = folder / "SKILL.md"
    if not manifest.is_file():
        return [f"{folder}: missing SKILL.md"]
    try:
        content = manifest.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        return [f"{manifest}: {exc}"]
    match = FRONTMATTER.match(content)
    if not match:
        return [f"{manifest}: missing or malformed YAML frontmatter"]
    try:
        metadata = yaml.load(match.group(1), Loader=UniqueKeyLoader)
    except (yaml.YAMLError, ValueError) as exc:
        return [f"{manifest}: invalid YAML: {exc}"]
    if not isinstance(metadata, dict):
        return [f"{manifest}: frontmatter must be a mapping"]
    errors = []
    name = metadata.get("name")
    if not isinstance(name, str) or not 1 <= len(name) <= 64 or not NAME.fullmatch(name):
        errors.append(f"{manifest}: name must be 1-64 lowercase letters/hyphens, without empty segments")
    if name != folder.name:
        errors.append(f"{manifest}: name must match parent directory")
    description = metadata.get("description")
    if not isinstance(description, str) or not 1 <= len(description.strip()) <= 1024:
        errors.append(f"{manifest}: description must be a nonempty string of at most 1024 characters")
    unknown = set(metadata) - STANDARD_FIELDS
    if unknown:
        errors.append(f"{manifest}: unsupported fields: {', '.join(sorted(unknown))}")
    for field in ("license", "compatibility", "allowed-tools"):
        if field in metadata and not isinstance(metadata[field], str):
            errors.append(f"{manifest}: {field} must be a string")
    if isinstance(metadata.get("compatibility"), str) and not 1 <= len(metadata["compatibility"]) <= 500:
        errors.append(f"{manifest}: compatibility must be 1-500 characters")
    if "metadata" in metadata and (
        not isinstance(metadata["metadata"], dict)
        or any(not isinstance(v, str) for v in metadata["metadata"].values())
    ):
        errors.append(f"{manifest}: metadata must map string keys to string values")
    if not content[match.end():].strip():
        errors.append(f"{manifest}: instructions are empty")
    return errors


def validate_links(document: Path, root: Path) -> list[str]:
    """Validate inline local file links, resolving relative to their document."""
    errors = []
    try:
        content = document.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        return [f"{document}: {exc}"]
    # Ignore illustrative links in fenced code blocks.
    content = re.sub(r"(?ms)^```[^\n]*\n.*?^```[^\n]*$", "", content)
    for match in INLINE_LINK.finditer(content):
        target = match.group(1).strip("<>")
        parsed = urlsplit(target)
        if parsed.scheme in ("http", "https", "mailto") or target.startswith("#"):
            continue
        if parsed.scheme or parsed.netloc or Path(parsed.path).is_absolute():
            errors.append(f"{document}: use a relative local link: {target}")
            continue
        resolved = (document.parent / unquote(parsed.path)).resolve()
        if not resolved.is_relative_to(root.resolve()):
            errors.append(f"{document}: local link leaves repository: {target}")
        elif not resolved.exists():
            errors.append(f"{document}: broken local link: {target}")
    return errors


def validate_repository(root: Path) -> tuple[list[str], int, int]:
    """Check the seven required skills, any additions, and guidance file links."""
    skills = root / ".agents" / "skills"
    folders = sorted(p for p in skills.iterdir() if p.is_dir()) if skills.is_dir() else []
    errors = [f"Missing required skill: {name}" for name in sorted(REQUIRED_SKILLS - {p.name for p in folders})]
    for folder in folders:
        errors.extend(validate_skill(folder))
    documents = list(skills.rglob("*.md")) if skills.is_dir() else []
    for relative in ("AGENTS.md", "README.md", "docs/AGENT_SKILLS_VALIDATION.md"):
        path = root / relative
        if path.is_file():
            documents.append(path)
        else:
            errors.append(f"Missing guidance: {relative}")
    for document in documents:
        errors.extend(validate_links(document, root))
    return errors, len(folders), len(documents)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[4], help="Repository root (defaults to this checkout)")
    args = parser.parse_args()
    errors, skills, documents = validate_repository(args.root.resolve())
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    if errors:
        print(f"Failed with {len(errors)} issue(s).", file=sys.stderr)
        return 1
    print(f"PASS: {skills} skills; metadata and inline local links in {documents} Markdown files.")
    print("Manual review still required: descriptions, API accuracy, scientific guidance, external links, and real activation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
