# Auto MPG Data Set

Data concerns city-cycle fuel consumption - revised from CMU StatLib
library.

## Usage

``` r
auto_mpg
```

## Format

A matrix containing 398 observations and 10 attributes.

- `mpg`:

  Miles per gallon of the engine. Predictor attribute

- `cylinders`:

  Number of cylinders in the engine

- `displacement`:

  Engine displacement

- `horsepower`:

  Horsepower of the car

- `weight`:

  Weight of the car (lbs)

- `acceleration`:

  Acceleration of the car (seconds taken for 0-60mph)

- `model_year`:

  Model year of the car in the 1900s

- `origin`:

  Car origin

- `make`:

  Car manufacturer

- `car_name`:

  Name of the car

## Source

http://archive.ics.uci.edu/ml/datasets/Auto+MPG

## References

Quinlan,R. (1993). Combining Instance-Based and Model-Based Learning. In
Proceedings on the Tenth International Conference of Machine Learning,
236-243, University of Massachusetts, Amherst. Morgan Kaufmann.

## Examples

``` r
data(auto_mpg)    # Lazy loading
```
