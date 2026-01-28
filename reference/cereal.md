# Breakfast Cereal Data Set

Nutritional information and manufacturer data for 70+ popular US
breakfast cereals

## Usage

``` r
cereal
```

## Format

A matrix containing 77 observations and 16 attributes.

- `name`:

  name of cereal

- `manuf`:

  manufacturer of cereal, coded into seven categories: "A" for American
  Home Food Products, "G" for General Mills, "K" for Kelloggs, "N" for
  Nabisco, "P" for Post, "Q" for Quaker Oats, and "R" for Ralston Purina

- `type`:

  cold or hot

- `calories`:

  calories per serving

- `protein`:

  grams of protein

- `fat`:

  grams of fat

- `sodium`:

  milligrams of sodium

- `fiber`:

  grams of dietary fiber

- `carbo`:

  grams of complex carbohydrates

- `sugars`:

  grams of sugars

- `potass`:

  milligrams of potassium

- `vitamins`:

  vitamins and minerals - 0, 25, or 100, indicating the typical
  percentage of FDA recommended

- `shelf`:

  display shelf (1, 2, or 3, counting from the floor)

- `weight`:

  weight in ounces of one serving

- `cups`:

  number of cups in one serving

- `rating`:

  a rating of the cereals

## Source

https://lib.stat.cmu.edu/datasets/1993.expo/

## References

Reza Mohammadi (2025). Data Science Foundations and Machine Learning
with R: From Data to Decisions.
<https://book-data-science-r.netlify.app>.

## Examples

``` r
data(cereal)
str(cereal)
#> 'data.frame':    77 obs. of  16 variables:
#>  $ name    : Factor w/ 77 levels "100% Bran","100% Natural Bran",..: 1 2 3 4 5 6 7 8 9 10 ...
#>  $ manuf   : Factor w/ 7 levels "A","G","K","N",..: 4 6 3 3 7 2 3 2 7 5 ...
#>  $ type    : Factor w/ 2 levels "cold","hot": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ calories: int  70 120 70 50 110 110 110 130 90 90 ...
#>  $ protein : int  4 3 4 4 2 2 2 3 2 3 ...
#>  $ fat     : int  1 5 1 0 2 2 0 2 1 0 ...
#>  $ sodium  : int  130 15 260 140 200 180 125 210 200 210 ...
#>  $ fiber   : num  10 2 9 14 1 1.5 1 2 4 5 ...
#>  $ carbo   : num  5 8 7 8 14 10.5 11 18 15 13 ...
#>  $ sugars  : int  6 8 5 0 8 10 14 8 6 5 ...
#>  $ potass  : int  280 135 320 330 -1 70 30 100 125 190 ...
#>  $ vitamins: int  25 0 25 25 25 25 25 25 25 25 ...
#>  $ shelf   : int  3 3 3 3 3 1 2 3 1 3 ...
#>  $ weight  : num  1 1 1 1 1 1 1 1.33 1 1 ...
#>  $ cups    : num  0.33 1 0.33 0.5 0.75 0.75 1 0.75 0.67 0.67 ...
#>  $ rating  : num  68.4 34 59.4 93.7 34.4 ...
```
