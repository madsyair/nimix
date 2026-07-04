# Clear the compiled-model cache

nimix reuses compiled NIMBLE models across fits that share an identical
model structure (see the `reuse` entry of `mcmcControl` in
[`nimixClust`](https://madsyair.github.io/nimix/reference/nimixClust.md)).
Compiled models are large; this empties the cache and releases them.

## Usage

``` r
nimixClearCache()
```

## Value

Invisibly, the number of cached compiled models that were removed.

## Examples

``` r
nimixClearCache()
```
