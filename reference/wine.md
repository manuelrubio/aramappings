# Wine Data Set

Chemical analysis to determine the origin of wines.

## Usage

``` r
wine
```

## Format

A matrix containing 178 observations and 14 attributes (including 1
classification attribute).

- `Class`:

  Class of cultivar

- `Alcohol`:

  Alcohol

- `Malic acid`:

  Malic acid

- `Ash`:

  Ash

- `Alcalinity of ash`:

  Alcalinity of ash

- `Magnesium`:

  Magnesium

- `Total phenols`:

  Total phenols

- `Flavanoids`:

  Flavanoids

- `Nonflavanoid phenols`:

  Nonflavanoid phenols

- `Proanthocyanins`:

  Proanthocyanins

- `Color intensity`:

  Color intensity

- `Hue`:

  Hue

- `OD280/OD315 of diluted wines`:

  OD280/OD315 of diluted wines

- `Proline`:

  Proline

## Source

http://archive.ics.uci.edu/dataset/109/wine

## References

Dua, D., Graff, C.: UCI Machine Learning Repository. University of
California, School of Information and Computer Science, Irvine, CA
(2019)

## Examples

``` r
data(wine)
X <- wine[,-1]
class <- wine[,1]
```
