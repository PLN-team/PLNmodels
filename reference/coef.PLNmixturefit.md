# Extract model coefficients

Extracts model coefficients from objects returned by
[`PLN()`](https://pln-team.github.io/PLNmodels/reference/PLN.md) and its
variants

## Usage

``` r
# S3 method for class 'PLNmixturefit'
coef(object, type = c("main", "means", "covariance", "mixture"), ...)
```

## Arguments

- object:

  an R6 object with class
  [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)

- type:

  type of parameter that should be extracted. Either "main" (default)
  for \$\$\Theta\$\$, "means" for \$\$\mu\$\$, "mixture" for \$\$\pi\$\$
  or "covariance" for \$\$\Sigma\$\$

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A matrix of coefficients extracted from the PLNfit model.

## See also

[`sigma.PLNmixturefit()`](https://pln-team.github.io/PLNmodels/reference/sigma.PLNmixturefit.md)

## Examples

``` r
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
myPLN <- PLNmixture(Abundance ~ 1 + offset(log(Offset)),
           data = trichoptera, control = PLNmixture_param(smoothing = "none"))  %>% getBestModel()
#> 
#>  Initialization...
#> 
#>  Adjusting 5 PLN mixture models.
#>  number of cluster = 1   number of cluster = 2   number of cluster = 3   number of cluster = 4   number of cluster = 5 
#>  Post-treatments
#>  DONE!
coef(myPLN) ## Theta - empty here
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#>      [,15] [,16] [,17]
coef(myPLN, type = "mixture") ## pi
#> [1] 0.04081633 0.16301704 0.63290133 0.16326531
coef(myPLN, type = "means") ## mu
#>       Intercept Intercept.1 Intercept.2 Intercept.3
#> Che -12.5029490 -12.8059740   -6.498698  -12.245351
#> Hyc  -8.5026525 -12.8059740   -7.655274   -5.654305
#> Hym  -6.2095447  -3.8209871   -2.400798   -3.226406
#> Hys -12.5029490 -12.8059740   -5.946180   -5.587155
#> Psy  -0.2098338  -0.2168007   -0.575118   -1.129249
#> Aga  -7.3759594  -3.4593157   -3.778450   -5.598457
#> Glo  -6.5250622 -12.8059740   -5.797673   -5.598457
#> Ath  -7.8104125  -7.0212752   -5.336241   -5.576372
#> Cea  -6.8782910 -12.8059740   -6.887291  -12.245351
#> Ced  -5.6488768  -4.7794263   -3.232171   -3.268467
#> Set  -3.6128155  -5.6111089   -4.050179   -3.119234
#> All  -8.5599770  -4.7457972   -4.161607   -5.612215
#> Han  -4.6791420  -4.8807234   -3.924914   -2.331123
#> Hfo  -3.5993025 -12.8059740   -6.043910   -3.581710
#> Hsp  -2.4358493  -5.1354675   -4.937336   -1.482940
#> Hve  -7.8104125  -7.6633804   -5.850092  -12.245351
#> Sta  -5.7700392  -2.8104490   -2.697348   -2.531135
coef(myPLN, type = "covariance") ## Sigma
#> [[1]]
#> 17 x 17 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 17 column names ‘Che’, ‘Hyc’, ‘Hym’ ... ]]
#>                                                                          
#> Che 0.2080685 .         .         .         .         .         .        
#> Hyc .         0.2080685 .         .         .         .         .        
#> Hym .         .         0.2080685 .         .         .         .        
#> Hys .         .         .         0.2080685 .         .         .        
#> Psy .         .         .         .         0.2080685 .         .        
#> Aga .         .         .         .         .         0.2080685 .        
#> Glo .         .         .         .         .         .         0.2080685
#> Ath .         .         .         .         .         .         .        
#> Cea .         .         .         .         .         .         .        
#> Ced .         .         .         .         .         .         .        
#> Set .         .         .         .         .         .         .        
#> All .         .         .         .         .         .         .        
#> Han .         .         .         .         .         .         .        
#> Hfo .         .         .         .         .         .         .        
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                                                          
#> Che .         .         .         .         .         .         .        
#> Hyc .         .         .         .         .         .         .        
#> Hym .         .         .         .         .         .         .        
#> Hys .         .         .         .         .         .         .        
#> Psy .         .         .         .         .         .         .        
#> Aga .         .         .         .         .         .         .        
#> Glo .         .         .         .         .         .         .        
#> Ath 0.2080685 .         .         .         .         .         .        
#> Cea .         0.2080685 .         .         .         .         .        
#> Ced .         .         0.2080685 .         .         .         .        
#> Set .         .         .         0.2080685 .         .         .        
#> All .         .         .         .         0.2080685 .         .        
#> Han .         .         .         .         .         0.2080685 .        
#> Hfo .         .         .         .         .         .         0.2080685
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                  
#> Che .         .         .        
#> Hyc .         .         .        
#> Hym .         .         .        
#> Hys .         .         .        
#> Psy .         .         .        
#> Aga .         .         .        
#> Glo .         .         .        
#> Ath .         .         .        
#> Cea .         .         .        
#> Ced .         .         .        
#> Set .         .         .        
#> All .         .         .        
#> Han .         .         .        
#> Hfo .         .         .        
#> Hsp 0.2080685 .         .        
#> Hve .         0.2080685 .        
#> Sta .         .         0.2080685
#> 
#> [[2]]
#> 17 x 17 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 17 column names ‘Che’, ‘Hyc’, ‘Hym’ ... ]]
#>                                                                          
#> Che 0.5198944 .         .         .         .         .         .        
#> Hyc .         0.5198944 .         .         .         .         .        
#> Hym .         .         0.5198944 .         .         .         .        
#> Hys .         .         .         0.5198944 .         .         .        
#> Psy .         .         .         .         0.5198944 .         .        
#> Aga .         .         .         .         .         0.5198944 .        
#> Glo .         .         .         .         .         .         0.5198944
#> Ath .         .         .         .         .         .         .        
#> Cea .         .         .         .         .         .         .        
#> Ced .         .         .         .         .         .         .        
#> Set .         .         .         .         .         .         .        
#> All .         .         .         .         .         .         .        
#> Han .         .         .         .         .         .         .        
#> Hfo .         .         .         .         .         .         .        
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                                                          
#> Che .         .         .         .         .         .         .        
#> Hyc .         .         .         .         .         .         .        
#> Hym .         .         .         .         .         .         .        
#> Hys .         .         .         .         .         .         .        
#> Psy .         .         .         .         .         .         .        
#> Aga .         .         .         .         .         .         .        
#> Glo .         .         .         .         .         .         .        
#> Ath 0.5198944 .         .         .         .         .         .        
#> Cea .         0.5198944 .         .         .         .         .        
#> Ced .         .         0.5198944 .         .         .         .        
#> Set .         .         .         0.5198944 .         .         .        
#> All .         .         .         .         0.5198944 .         .        
#> Han .         .         .         .         .         0.5198944 .        
#> Hfo .         .         .         .         .         .         0.5198944
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                  
#> Che .         .         .        
#> Hyc .         .         .        
#> Hym .         .         .        
#> Hys .         .         .        
#> Psy .         .         .        
#> Aga .         .         .        
#> Glo .         .         .        
#> Ath .         .         .        
#> Cea .         .         .        
#> Ced .         .         .        
#> Set .         .         .        
#> All .         .         .        
#> Han .         .         .        
#> Hfo .         .         .        
#> Hsp 0.5198944 .         .        
#> Hve .         0.5198944 .        
#> Sta .         .         0.5198944
#> 
#> [[3]]
#> 17 x 17 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 17 column names ‘Che’, ‘Hyc’, ‘Hym’ ... ]]
#>                                                                          
#> Che 0.6453444 .         .         .         .         .         .        
#> Hyc .         0.6453444 .         .         .         .         .        
#> Hym .         .         0.6453444 .         .         .         .        
#> Hys .         .         .         0.6453444 .         .         .        
#> Psy .         .         .         .         0.6453444 .         .        
#> Aga .         .         .         .         .         0.6453444 .        
#> Glo .         .         .         .         .         .         0.6453444
#> Ath .         .         .         .         .         .         .        
#> Cea .         .         .         .         .         .         .        
#> Ced .         .         .         .         .         .         .        
#> Set .         .         .         .         .         .         .        
#> All .         .         .         .         .         .         .        
#> Han .         .         .         .         .         .         .        
#> Hfo .         .         .         .         .         .         .        
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                                                          
#> Che .         .         .         .         .         .         .        
#> Hyc .         .         .         .         .         .         .        
#> Hym .         .         .         .         .         .         .        
#> Hys .         .         .         .         .         .         .        
#> Psy .         .         .         .         .         .         .        
#> Aga .         .         .         .         .         .         .        
#> Glo .         .         .         .         .         .         .        
#> Ath 0.6453444 .         .         .         .         .         .        
#> Cea .         0.6453444 .         .         .         .         .        
#> Ced .         .         0.6453444 .         .         .         .        
#> Set .         .         .         0.6453444 .         .         .        
#> All .         .         .         .         0.6453444 .         .        
#> Han .         .         .         .         .         0.6453444 .        
#> Hfo .         .         .         .         .         .         0.6453444
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                  
#> Che .         .         .        
#> Hyc .         .         .        
#> Hym .         .         .        
#> Hys .         .         .        
#> Psy .         .         .        
#> Aga .         .         .        
#> Glo .         .         .        
#> Ath .         .         .        
#> Cea .         .         .        
#> Ced .         .         .        
#> Set .         .         .        
#> All .         .         .        
#> Han .         .         .        
#> Hfo .         .         .        
#> Hsp 0.6453444 .         .        
#> Hve .         0.6453444 .        
#> Sta .         .         0.6453444
#> 
#> [[4]]
#> 17 x 17 sparse Matrix of class "dgCMatrix"
#>   [[ suppressing 17 column names ‘Che’, ‘Hyc’, ‘Hym’ ... ]]
#>                                                                          
#> Che 0.2412511 .         .         .         .         .         .        
#> Hyc .         0.2412511 .         .         .         .         .        
#> Hym .         .         0.2412511 .         .         .         .        
#> Hys .         .         .         0.2412511 .         .         .        
#> Psy .         .         .         .         0.2412511 .         .        
#> Aga .         .         .         .         .         0.2412511 .        
#> Glo .         .         .         .         .         .         0.2412511
#> Ath .         .         .         .         .         .         .        
#> Cea .         .         .         .         .         .         .        
#> Ced .         .         .         .         .         .         .        
#> Set .         .         .         .         .         .         .        
#> All .         .         .         .         .         .         .        
#> Han .         .         .         .         .         .         .        
#> Hfo .         .         .         .         .         .         .        
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                                                          
#> Che .         .         .         .         .         .         .        
#> Hyc .         .         .         .         .         .         .        
#> Hym .         .         .         .         .         .         .        
#> Hys .         .         .         .         .         .         .        
#> Psy .         .         .         .         .         .         .        
#> Aga .         .         .         .         .         .         .        
#> Glo .         .         .         .         .         .         .        
#> Ath 0.2412511 .         .         .         .         .         .        
#> Cea .         0.2412511 .         .         .         .         .        
#> Ced .         .         0.2412511 .         .         .         .        
#> Set .         .         .         0.2412511 .         .         .        
#> All .         .         .         .         0.2412511 .         .        
#> Han .         .         .         .         .         0.2412511 .        
#> Hfo .         .         .         .         .         .         0.2412511
#> Hsp .         .         .         .         .         .         .        
#> Hve .         .         .         .         .         .         .        
#> Sta .         .         .         .         .         .         .        
#>                                  
#> Che .         .         .        
#> Hyc .         .         .        
#> Hym .         .         .        
#> Hys .         .         .        
#> Psy .         .         .        
#> Aga .         .         .        
#> Glo .         .         .        
#> Ath .         .         .        
#> Cea .         .         .        
#> Ced .         .         .        
#> Set .         .         .        
#> All .         .         .        
#> Han .         .         .        
#> Hfo .         .         .        
#> Hsp 0.2412511 .         .        
#> Hve .         0.2412511 .        
#> Sta .         .         0.2412511
#> 
```
