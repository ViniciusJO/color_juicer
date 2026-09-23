# Color juicer [WIP]

A small, opinionated command-line tool that extracts a color palette from an image using [k-means](https://en.wikipedia.org/wiki/K-means_clustering) clustering, then derives a contrasting/complementary color for every extracted color. It can also render those colors into arbitrary template files (theme files for a status bar, a terminal, a window manager, etc.), mirroring a whole directory tree in the process.

Built in [Zig](https://ziglang.org/) (0.16), using [stb_image/stb_image_resize2/stb_image_write](https://github.com/nothings/stb) for image decoding/downscaling and [zig-args](https://github.com/ikskuh/zig-args) for CLI parsing.

> The result is **not deterministic**: the initial centroids are chosen with a randomized k-means++ seeding step, so re-running the tool on the same image can produce a (usually similar, not identical) palette every time.

## How it works

1. **Load & downsample** — the input image is decoded with `stb_image` (forced to 3 channels, RGB) and downsampled by a fixed factor (`5x` on each axis) with `stb_image_resize2`. Downsampling both speeds up clustering and smooths out noise/compression artifacts.
2. **k-means++ initialization** — the first centroid is picked at random from the downsampled pixels; each subsequent centroid is picked with probability proportional to its squared distance from the closest already-chosen centroid. This spreads the initial centroids out instead of risking several of them landing on near-identical colors.
3. **Lloyd's algorithm** — pixels are repeatedly repartitioned to their nearest centroid (Euclidean distance in RGB space) and centroids are recomputed as the average color of their partition, until no centroid moves more than a small threshold (`0.01`).
4. **Complementary colors** — every extracted color is converted to [OKLCH](https://bottosson.github.io/posts/oklab/) and rotated 180° in hue to get a perceptually-reasonable complementary color.
5. **Contrast-aware sorting** — the palette is sorted by how strongly each extracted color contrasts against a reference color (`--contrast`, `#000000` by default): colors are ranked by a combination of raw RGB distance and the delta in OKLCH lightness/chroma against the reference. This is what makes the *first* colors in the output the ones that "pop" the most against your reference (e.g. your terminal/bar background).
6. **Output** — the sorted palette (primaries `p1..pN` and their complementaries `c1..cN`) is printed to stdout, optionally written to `$HOME/.cache/juiced.color`, and optionally rendered into a tree of template files.

## Building

Requires **Zig 0.16**. The repo pins this via [`mise`](https://mise.jdx.dev/) (`mise.toml` -> `zig = "0.16"`); if you're not using `mise`, just make sure `zig version` reports `0.16.x`.

```sh
zig build                       # debug build -> zig-out/bin/color_juicer
zig build -Doptimize=ReleaseFast
zig build run -- <args...>      # build + run in one step
```

`build.zig.zon` pulls in the `zig-args` dependency automatically; the `stb_*` headers are vendored under `src/stb/` and compiled/translated as part of the build (no system image libraries required, everything links statically).

CI (`.github/workflows/zigbuild.yml`) builds a `ReleaseFast` static binary on every push/PR to the `release` branch and publishes it as a GitHub release artifact.

## Usage

```
color_juicer <image-path> [options]
```

```
Options:
  -#, --colors_count <n>         number of colors to output (defaults to 6)
  -c, --contrast <hex>           color on which to calculate contrast (defaults to #000000)
  -t, --templates_path <path>    path for the template folder
  -o, --output_reference_path <path>
                                 reference base path to resolve templates (defaults to $HOME)
  -n, --no_cache_output          prevent generation of the default $HOME/.cache/juiced.color file
  -h, --help                     displays help message
```

Notes on argument syntax (inherited from `zig-args`):

- Short options that take a value need a space (or `=`) between the flag and the value: `-# 4` or `-#=4`; something like `-#4` will *not* parse.
- Long options accept either a space or `=`: `--colors_count 4` or `--colors_count=4`.
- `--contrast` accepts `#RGB`, `#RGBA`, `#RRGGBB` or `#RRGGBBAA` hex strings (case-insensitive). Anything else is rejected with an error and a non-zero exit code.

### Examples

Extract the default 6 colors, contrasted against black, write `$HOME/.cache/juiced.color`:

```sh
color_juicer photo.png
```

Extract 4 colors, contrasted against white, without touching the cache file:

```sh
color_juicer photo.png -# 4 -c "#ffffff" -n
```

Render a set of templates using the extracted palette:

```sh
color_juicer photo.png -t ~/.config/color_juicer/template -o ~/.cache/color_juicer
```

## Output

Using the reference image [test.png](test.png):

![test.png](test.png)

### Console

For every run, `color_juicer` prints the raw list of extracted colors (in clustering/partition-size order) followed by the contrast reference color and the final, contrast-sorted palette:


> Image_path: test.png
> 
> Input: Image{ .width = 3440, .height = 1440, .channels = 3 }
> Down:  Image{ .width = 688, .height = 288, .channels = 3 }
> 
> Extracted_colors:
> 
>  ![#222D36](https://placehold.co/15x15/222d36/222d36.png) #222D36
>  ![#EFD6A4](https://placehold.co/15x15/efd6a4/efd6a4.png) #EFD6A4
>  ![#AAB4A1](https://placehold.co/15x15/aab4a1/aab4a1.png) #AAB4A1
>  ![#427D95](https://placehold.co/15x15/427d95/427d95.png) #427D95
>  ![#5C5F4C](https://placehold.co/15x15/5c5f4c/5c5f4c.png) #5C5F4C
>  ![#A28762](https://placehold.co/15x15/a28762/a28762.png) #A28762
> 
> Contrast_color: ![#000000](https://placehold.co/15x15/000000/000000.png) #000000
> 
> Pallete:
> 
> p1: ![#EFD6A4](https://placehold.co/15x15/efd6a4/efd6a4.png) #EFD6A4
> p2: ![#A28762](https://placehold.co/15x15/a28762/a28762.png) #A28762
> p3: ![#AAB4A1](https://placehold.co/15x15/aab4a1/aab4a1.png) #AAB4A1
> p4: ![#427D95](https://placehold.co/15x15/427d95/427d95.png) #427D95
> p5: ![#5C5F4C](https://placehold.co/15x15/5c5f4c/5c5f4c.png) #5C5F4C
> p6: ![#222D36](https://placehold.co/15x15/222d36/222d36.png) #222D36
> 
> c1: ![#C2D9FF](https://placehold.co/15x15/c2d9ff/c2d9ff.png) #C2D9FF
> c2: ![#738EB1](https://placehold.co/15x15/738eb1/738eb1.png) #738EB1
> c3: ![#B6ABBE](https://placehold.co/15x15/b6abbe/b6abbe.png) #B6ABBE
> c4: ![#986750](https://placehold.co/15x15/986750/986750.png) #986750
> c5: ![#5E5A6C](https://placehold.co/15x15/5e5a6c/5e5a6c.png) #5E5A6C
> c6: ![#342920](https://placehold.co/15x15/342920/342920.png) #342920
> 
> juiced_colors: /home/user/.cache/juiced.color

(colors are printed with real ANSI 24-bit color escapes in an actual terminal — the swatches above are just their hex values for reference.)

`pN` are the extracted palette colors, sorted from most- to least-contrasting against `--contrast`; `cN` are each `pN`'s OKLCH complementary color, at the same index.

### `$HOME/.cache/juiced.color`

Unless `-n`/`--no_cache_output` is passed (and `$HOME` is set), the same palette is written to `$HOME/.cache/juiced.color` in a simple, greppable, line-oriented format — one color per line, primaries first, then complementaries:

```
p1 #EFD6A4
p2 #A28762
p3 #AAB4A1
p4 #427D95
p5 #5C5F4C
p6 #222D36
c1 #C2D9FF
c2 #738EB1
c3 #B6ABBE
c4 #986750
c5 #5E5A6C
c6 #342920
```

This file is meant to be sourced/parsed by other scripts (shell, polybar, i3/sway config, etc.) that don't need the full template engine below.

## Template engine

Pass `-t/--templates_path <dir>` to render an entire directory tree of templates with the extracted palette. Each file under `<dir>` (recursively) is copied to the same relative path under `-o/--output_reference_path` (defaults to `$HOME`), with every `%%pN%%` / `%%cN%%` placeholder replaced by the corresponding color's `#RRGGBB` hex string. Directories that don't exist yet on the output side are created; files with the same relative path are overwritten.

Placeholder syntax: `%%` followed by `p` or `c` and a 1-based index — e.g. `%%p1` is the most-contrasting primary color, `%%c2` is the complementary of the second primary color. `N` must be within `1..colors_count`; out-of-range or malformed placeholders are silently replaced with an empty string.

Example — given `~/tpl/colors.conf`:

```ini
[colors]
background = %%p1
foreground = %%c1
accent     = %%p2
```

running:

```sh
color_juicer photo.png -t ~/tpl -o ~/.cache/color_juicer
```

produces `~/.cache/color_juicer/colors.conf` (directory created automatically) with the placeholders substituted, e.g.:

```ini
[colors]
background = #EFD6A4
foreground = #C2D9FF
accent     = #A28762
```

While generating, the tool prints the resolved template/output paths and a small colored tree of every file/directory it visits, so you can confirm what was written before trusting the result.

## Project layout

```
src/
  main.zig               CLI entry point: argument parsing, image loading/
                          downsampling, k-means clustering, sorting, output
  color.zig               RGBA <-> OKLCH color conversions, complementary
                          color, hex formatting/parsing helpers
  vec.zig                 small fixed-size color vectors (Vec3/Vec4) and
                          distance helpers, plus safe byte-buffer -> Vec3
                          reinterpretation (bytes_as_vec3)
  pallete.zig              the extracted palette (primaries/complementaries)
                          and its console/file printers
  parser.zig               the %%pN%%/%%cN%% template engine and recursive
                          directory-tree renderer
  usage_description.txt    long-form help text embedded into the --help output
  stb/                     vendored stb_image / stb_image_resize2 /
                          stb_image_write headers, plus hand-written extern
                          bindings (stb/bindings.zig) as a lighter
                          alternative to full @cImport translation
  .ignore/                 local scratch templates/output used during
                          development (git-ignored, not part of the tool)
```

## Known limitations

- Non-deterministic output (see above) — re-run if you don't like the result.
- Only image formats supported by `stb_image` are readable (PNG, JPEG, BMP, GIF, PSD, TGA, HDR, PIC, PNM — see the [stb_image README](https://github.com/nothings/stb/blob/master/stb_image.h) for the exact list/caveats).
- The template placeholder scanner does not flush a placeholder that is still open at end-of-file (a `%%p1` with no closing `%%` before EOF is dropped instead of emitted literally).
- Short options with a value (`-#`, `-c`, `-t`, `-o`) must be given their value as a separate token (or with `=`); a glued form like `-#4` fails to parse.

## TBD

- [ ] To reduce binary size, try using binding only for used STB functions
- [ ] Correctly contrast calculation and sorting
- [ ] More output/template presets out of the box (rofi, dunst, GTK, Qt, ...).
- [ ] General code cleanup of `main.zig` (still carries a fair amount of commented-out exploration code from earlier iterations).
