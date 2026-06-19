# Extract variance-covariance of residuals 'Sigma'

Extract the variance-covariance matrix of the residuals, usually noted
\$\$\Sigma\$\$ in PLN models. This captures the correlation between the
species in the latent space. or PLNmixture, it is a weighted mean of the
variance-covariance matrices of each component.

## Usage

``` r
# S3 method for class 'PLNmixturefit'
sigma(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`PLNmixturefit`](https://pln-team.github.io/PLNmodels/reference/PLNmixturefit.md)

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A semi definite positive matrix of size p, assuming there are p species
in the model.

## See also

[`coef.PLNmixturefit()`](https://pln-team.github.io/PLNmodels/reference/coef.PLNmixturefit.md)
for other ways to access \$\$\Sigma\$\$.

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
sigma(myPLN) ## Sigma
#> [[1]]
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
#> [[2]]
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
#> [[3]]
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
#> [[4]]
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
```
