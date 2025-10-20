import subprocess
import tempfile
import os
import argparse
from generate_tiff_offsets import generate_tiff_offsets

def bf_to_jp2k(input_path: str, output_path: str, lossless: bool = False, quality: int = 10):
    """
    Convert a Bioformats-compatible file to a JP2K OME-TIFF file.
    Requires bioformats2raw & raw2ometiff to be installed.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # use bioformats2raw to convert the file to a raw file
        raw_path = os.path.join(tmpdir, "raw.raw")
        subprocess.run(["bioformats2raw", input_path, raw_path], check=True)
        # use raw2ometiff to convert the raw file to a JP2K OME-TIFF file
        subprocess.run([
            "raw2ometiff", raw_path, output_path,
            "--compression", "JPEG-2000" if lossless else "JPEG-2000 Lossy",
            "--quality", str(quality),
        ], check=True)
        generate_tiff_offsets(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", type=str)
    parser.add_argument("output_path", type=str)
    parser.add_argument("--lossless", action="store_true")
    parser.add_argument("--quality", type=int, default=10)
    args = parser.parse_args()
    print(f"Converting {args.input_path} to JP2K OME-TIFF...")
    bf_to_jp2k(args.input_path, args.output_path, args.lossless, args.quality)
    print(f"Converted {args.input_path} to {args.output_path}")