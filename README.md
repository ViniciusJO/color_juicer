# Color juicer

Uses the kmeans algorithm to extract colors from image.

Made using [stb_image](https://github.com/nothings/stb) utilities;

It is extremly opinionated to get the most brigth colors to contrast with the black background of my polybar. The code also calculates the color that most contrast with the ones extracted in order to use as accent in the system.

Configured to generate the following files:

- `$HOME/.cache/dyn_colors.ini`: to use in polybar (exposes colors pprim, psec, pterc, pcont);
- `$HOME/.cache/i3_colors`: to configure the borders of focused containers in i3 with the color that most contrast with the extracted ones

## TBD:

Generate files to stylize:

- rofi
- dunst
- gtk
- qt
