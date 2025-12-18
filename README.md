# Color juicer

Uses the [k-means](https://pt.wikipedia.org/wiki/K-means) algorithm to extract colors from image.

Made using [stb_image](https://github.com/nothings/stb) utilities;

It is extremly opinionated to get the most brigth colors to contrast with the black background of my polybar. The code also calculates the color that most contrast with the ones extracted, in order to use as accent in the system.

\# The result is not deterministic !!! 

Configured to generate the following files:

- `$HOME/.cache/dyn_colors.ini`: to use in [polybar](https://github.com/polybar/polybar) (exposes colors pprim, psec, pterc, pcont);
- `$HOME/.cache/i3_colors`: to configure the borders of focused containers in [i3](https://github.com/i3/i3) with the color that most contrast with the extracted ones

## Output example

For the image:

![test.png](test.png)

extracts:

- pprim:  ![#ecd5a5](https://placehold.co/15x15/ecd5a5/ecd5a5.png) `#ecd5a5`
- psec:   ![#427b92](https://placehold.co/15x15/427b92/427b92.png) `#427b92`
- pterc:  ![#ada58c](https://placehold.co/15x15/ada58c/ada58c.png) `#ada58c`
- pcont:  ![#716365](https://placehold.co/15x15/716365/716365.png) `#716365`
 
and generates:
 
- on console

```
pprim: #ECD5A5
psec:  #427B92
pterc: #ADA58C
pcont: #716365
```

- on `$HOME/.cache/dyn_colors.ini`: 

```
[dyn_colors]
pprim = #ECD5A5
psec = #427B92
pterc = #ADA58C
pcont = #716365
cprim = #C7A795
csec = #5A6387
```

- on `$HOME/.cache/i3_colors`: 

```
client.focused #716365 #716365 #000000 #ECD5A5 #716365
```


## TBD:

- Cleanup code
- Generate files to stylize:
    - rofi
    - dunst
    - gtk
    - qt
