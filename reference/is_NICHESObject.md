# Test Whether an Object Is a NICHESObject

Convenience wrapper around `inherits(x, "NICHESObject")`. Useful in
downstream pipelines that need to validate inputs before calling
slot-aware helpers.

## Usage

``` r
is_NICHESObject(x)
```

## Arguments

- x:

  Any R object.

## Value

Logical scalar: `TRUE` if `x` inherits from class `"NICHESObject"`,
otherwise `FALSE`.

## Examples

``` r
is_NICHESObject(list())   # FALSE
#> [1] FALSE
is_NICHESObject(1:10)     # FALSE
#> [1] FALSE
```
