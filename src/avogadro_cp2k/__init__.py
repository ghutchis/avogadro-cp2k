#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

"""Avogadro plugin that reads CP2K main output files into Chemical JSON."""

import argparse
import json
import sys
import traceback

from .utils import decode, decompress


def main(argv=None):
    # Avogadro calls the plugin as:
    #   avogadro-cp2k <identifier> [--read] [--write] [--lang <locale>] [--debug]
    # Two input modes exist and both are accepted:
    #   - "file" mode (what pyproject.toml asks for): stdin carries a small JSON
    #     object such as {"operation": "read", "filename": "/path/to/run.out"}.
    #   - "string" mode: stdin carries the raw bytes of the file itself.
    # Either way stdout must be nothing but the bare CJSON document.
    parser = argparse.ArgumentParser()
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--read", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)

    # Read stdin as bytes so the platform locale never decides the encoding.
    raw = sys.stdin.buffer.read()

    try:
        filename = None
        payload = None
        # A CP2K output never starts with "{"; only then look for the envelope.
        if raw.lstrip()[:1] == b"{":
            try:
                payload = json.loads(decode(raw))
            except ValueError:
                payload = None
        if isinstance(payload, dict) and ("filename" in payload or "file" in payload):
            filename = payload.get("filename", payload.get("file"))
            with open(filename, "rb") as handle:
                raw_content = handle.read()
        else:
            raw_content = raw

        content = decode(decompress(raw_content))

        if args.debug:
            print(
                f"avogadro-cp2k: feature={args.feature} filename={filename!r} chars={len(content)}",
                file=sys.stderr,
            )

        output = None
        if args.feature == "cp2k" and args.read:
            from .reader import read_cp2k

            cjson = read_cp2k(content, filename=filename, debug=args.debug)
            output = json.dumps(cjson, ensure_ascii=True, allow_nan=False)
    except Exception as exc:  # noqa: BLE001 - report cleanly, never crash
        print(f"avogadro-cp2k: error reading input: {exc}", file=sys.stderr)
        if args.debug:
            traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    if output is None:
        print(
            f"avogadro-cp2k: no output produced for feature={args.feature!r}",
            file=sys.stderr,
        )
        sys.exit(1)

    sys.stdout.buffer.write(output.encode("utf-8"))
    sys.stdout.buffer.flush()
