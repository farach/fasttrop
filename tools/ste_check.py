#!/usr/bin/env python3
"""Find possible ASD-STE100 problems in the fasttrop documentation.

The script examines the rules of ASD-STE100 Issue 9 that a program can count.
It does not contain the STE dictionary. Thus, it cannot find all problems, and
some of its findings can be correct STE. A person who knows the standard must
examine each finding.

Errors:
- Rule 8.1: a semicolon.
- Rules 5.1 and 6.3: a sentence that is too long. The limit is 20 words for an
  instruction and 25 words for a description.
- Rule 6.6: a paragraph with more than six sentences.
- Rule 4.2: a contraction.
- Rule 3.2: a perfect tense or a progressive tense.
- Rule 1.1: a word from a list of words that STE does not approve.

Items to examine:
- Rule 3.5: a word that ends in "-ing" and that is not in the list of
  permitted words.
- Rule 3.6: a possible passive construction.
- Rule 9.3: a possible phrasal verb.
- Rule 1.3: a word that STE approves only with some meanings.

Usage, in the package folder:
    python tools/ste_check.py            # the package documentation
    python tools/ste_check.py FILE ...   # other files
The exit status is 1 if the script finds an error.
"""
import re
import sys
from pathlib import Path

DEFAULT_FILES = [
    "README.md", "NEWS.md", "cran-comments.md", "DESCRIPTION",
    "inst/bench/README.md", "inst/validation/README.md",
    "R/*.R", "src/trop_core.h",
]

# Words that end in "-ing" and that are not verb forms, or that are
# technical nouns in this package.
ING_OK = {
    "during", "nothing", "something", "anything", "everything", "thing",
    "things", "string", "strings", "bring", "king", "ring", "spring", "tuning",
    "throttling", "programming",
}
NOT_APPROVED = {
    "ensure", "ensures", "ensured", "utilize", "utilizes", "utilized",
    "utilise", "require", "requires", "required", "requirement",
    "requirements", "provide", "provides", "provided", "perform", "performs",
    "performed", "implement", "implements", "implemented", "facilitate",
    "leverage", "via", "obtain", "obtains", "obtained", "additional",
    "additionally", "sufficient", "numerous", "commence", "terminate",
    "assist", "indicate", "indicates", "verify", "verifies", "verified",
    "determine", "determines", "identical", "various", "prior", "regarding",
    "respectively", "should", "may", "might", "could", "would", "shall",
    "e.g.", "i.e.", "etc.", "furthermore", "moreover", "however",
}
REVIEW_WORDS = {
    "about": "means 'concerning', not 'approximately'",
    "since": "only for time, not for 'because'",
    "while": "only for time",
    "within": "use 'in' if possible",
    "per": "use 'for each' or 'in each'",
    "allow": "use 'let' or 'can'", "allows": "use 'let' or 'can'",
    "check": "noun only, or 'make sure', 'examine'",
    "following": "use 'these' or 'that follow'",
    "need": "use 'must' or 'necessary'", "needs": "use 'must' or 'necessary'",
}
PHRASAL = [
    "set up", "carry out", "find out", "point out", "look at", "look up",
    "turn off", "turn on", "back up", "give back", "hand back", "opt out",
    "make up", "figure out", "fill in", "rule out", "end up", "come up",
    "shut down", "speed up", "slow down", "put back", "run out", "pick up",
    "break down", "build up", "deal with", "take place", "go back",
]
IMPERATIVES = {
    "install", "run", "use", "set", "make", "do", "select", "read", "open",
    "type", "add", "remove", "write", "put", "keep", "get", "start", "stop",
    "change", "compare", "examine", "download", "load", "save", "delete",
    "update", "refer", "give", "find", "call", "copy", "clone",
}
ABBREV = ["Eq.", "Eqs.", "Fig.", "No.", "et al.", "vs.", "St."]
SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9`\"'(\[])")


def md_blocks(text):
    """Give the prose paragraphs of a Markdown file."""
    text = re.sub(r"<!--.*?-->", " ", text, flags=re.S)
    text = re.sub(r"```.*?```", "\n\n", text, flags=re.S)
    blocks, cur = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("|") or s.startswith("[!["):
            if cur:
                blocks.append(" ".join(cur))
                cur = []
            continue
        if re.match(r"^([-*+]|\d+\.)\s+", s):
            if cur:
                blocks.append(" ".join(cur))
            cur = [re.sub(r"^([-*+]|\d+\.)\s+", "", s)]
            continue
        cur.append(s)
    if cur:
        blocks.append(" ".join(cur))
    return blocks


def roxygen_blocks(text):
    """Give the prose paragraphs of the roxygen comments in an R file."""
    blocks, cur, skip = [], [], False
    first_line_of_block = True
    for line in text.splitlines():
        if not line.startswith("#'"):
            if cur:
                blocks.append(" ".join(cur))
            cur, skip, first_line_of_block = [], False, True
            continue
        s = line[2:].strip()
        if first_line_of_block:
            first_line_of_block = False  # The title is a heading.
            continue
        if not s:
            if cur:
                blocks.append(" ".join(cur))
            cur = []
            continue
        tag = re.match(r"^@(\w+)\s*(.*)$", s)
        if tag:
            if cur:
                blocks.append(" ".join(cur))
            cur = []
            name, rest = tag.group(1), tag.group(2)
            skip = name in {"examples", "export", "noRd", "keywords", "useDynLib",
                            "rdname", "name", "inheritParams", "importFrom",
                            "references", "seealso", "family"}
            if skip:
                continue
            if name == "param":
                rest = re.sub(r"^\S+\s*", "", rest)
            if name == "section":
                rest = re.sub(r"^[^:]*:\s*", "", rest)
            s = rest
            if not s:
                continue
        if skip:
            continue
        if re.match(r"^[-*]\s+", s):
            if cur:
                blocks.append(" ".join(cur))
            cur = [re.sub(r"^[-*]\s+", "", s)]
            continue
        cur.append(s)
    if cur:
        blocks.append(" ".join(cur))
    return blocks


def description_blocks(text):
    m = re.search(r"^Description:(.*?)(?=^\S)", text, flags=re.S | re.M)
    return [" ".join(m.group(1).split())] if m else []


def c_comment_blocks(text):
    blocks = []
    for m in re.finditer(r"/\*(.*?)\*/", text, flags=re.S):
        lines = [re.sub(r"^\s*\*\s?", "", ln).strip() for ln in m.group(1).splitlines()]
        para = []
        for ln in lines + [""]:
            if ln:
                para.append(ln)
            elif para:
                blocks.append(" ".join(para))
                para = []
    return blocks


def clean(sentence_block):
    s = re.sub(r"`[^`]*`", " CODE ", sentence_block)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)
    s = re.sub(r"https?://\S+|<doi:[^>]*>|\\doi\{[^}]*\}", " URL ", s)
    s = re.sub(r"\*\*|__|\*", "", s)
    s = re.sub(r"\[([^\]]+)\]", r"\1", s)
    return s


def sentences(block):
    s = clean(block)
    for i, a in enumerate(ABBREV):
        s = s.replace(a, a.replace(".", "<DOT%d>" % i))
    parts = [p.strip() for p in SENT_SPLIT.split(s) if p.strip()]
    out = []
    for p in parts:
        for i, a in enumerate(ABBREV):
            p = p.replace(a.replace(".", "<DOT%d>" % i), a)
        out.append(p)
    return out


def word_count(sentence):
    s = re.sub(r"\([^()]*\)", " PAREN ", sentence)  # Rule 8.5
    s = re.sub(r'"[^"]*"', " QUOTE ", s)              # Rule 8.6
    return len([w for w in s.split() if re.search(r"[A-Za-z0-9]", w)])


def is_instruction(sentence):
    words = re.findall(r"[A-Za-z']+", sentence)
    if not words:
        return False
    if words[0].lower() in IMPERATIVES:
        return True
    if words[0].lower() in {"if", "when", "to", "before", "after", "on", "for", "with"}:
        tail = sentence.split(",", 1)
        if len(tail) == 2:
            nxt = re.findall(r"[A-Za-z']+", tail[1])
            return bool(nxt) and nxt[0].lower() in IMPERATIVES
    return False


def check_block(block, where, errors, reviews):
    sents = sentences(block)
    if len(sents) > 6:
        errors.append((where, "6.6", "paragraph has %d sentences" % len(sents), block[:90]))
    for sent in sents:
        low = " " + sent.lower() + " "
        if ";" in sent:
            errors.append((where, "8.1", "semicolon", sent))
        n = word_count(sent)
        limit = 20 if is_instruction(sent) else 25
        if n > limit:
            errors.append((where, "5.1/6.3", "%d words (limit %d)" % (n, limit), sent))
        if re.search(r"\b\w+n't\b|\b\w+'(re|ve|ll|d)\b|\b(it|that|there|here|what|let)'s\b", sent, re.I):
            errors.append((where, "4.2", "contraction", sent))
        if re.search(r"\b(has|have|had)\s+(been|\w+ed|done|made|given|shown|written|taken|found|seen|run)\b", low):
            errors.append((where, "3.2", "perfect tense", sent))
        if re.search(r"\b(is|are|was|were|be|been)\s+\w+ing\b", low) and not re.search(r"\b(is|are|was|were|be|been)\s+(during|nothing|something|tuning)\b", low):
            errors.append((where, "3.2", "progressive tense", sent))
        for w in re.findall(r"[A-Za-z][A-Za-z.]*[A-Za-z.]|[A-Za-z]", sent):
            wl = w.lower().rstrip(".")
            if w.lower() in NOT_APPROVED or wl in NOT_APPROVED:
                errors.append((where, "1.1", "not approved: %s" % w, sent))
            if wl in REVIEW_WORDS:
                reviews.append((where, "1.3", "%s: %s" % (wl, REVIEW_WORDS[wl]), sent))
            if re.fullmatch(r"[a-z]+ings?", wl) and wl.rstrip("s") not in ING_OK and wl not in ING_OK:
                reviews.append((where, "3.5", "-ing word: %s" % w, sent))
        if re.search(r"\b(is|are|was|were|be|been|being)\s+(\w+ed|given|shown|written|done|made|taken|found|known|run|set|put|seen|built|chosen|drawn|held|kept)\b", low):
            reviews.append((where, "3.6", "possible passive", sent))
        for pv in PHRASAL:
            verb, part = pv.split(" ", 1)
            if re.search(r"\b%s\w*\s+%s\b" % (re.escape(verb[:-1] if verb.endswith("e") else verb), re.escape(part)), low):
                reviews.append((where, "9.3", "phrasal verb: %s" % pv, sent))


def blocks_for(path):
    text = Path(path).read_text(encoding="utf-8")
    name = Path(path).name
    if name.endswith(".md"):
        return md_blocks(text)
    if name == "DESCRIPTION":
        return description_blocks(text)
    if name.endswith(".R"):
        return roxygen_blocks(text)
    if name.endswith((".h", ".c")):
        return c_comment_blocks(text)
    return md_blocks(text)


def main(argv):
    patterns = argv or DEFAULT_FILES
    files = []
    for p in patterns:
        matches = sorted(Path(".").glob(p)) if any(c in p for c in "*?[") else [Path(p)]
        files.extend(m for m in matches if m.exists())
    errors, reviews, n_sent = [], [], 0
    for f in files:
        for k, block in enumerate(blocks_for(f)):
            n_sent += len(sentences(block))
            check_block(block, "%s, paragraph %d" % (f.as_posix(), k + 1), errors, reviews)
    for kind, items in (("ERROR", errors), ("EXAMINE", reviews)):
        for where, rule, msg, text in items:
            print("%s [rule %s] %s: %s\n    %s" % (kind, rule, where, msg, text[:160]))
    print("\n%d files, %d sentences, %d errors, %d items to examine."
          % (len(files), n_sent, len(errors), len(reviews)))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
