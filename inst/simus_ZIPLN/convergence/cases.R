## Real data cases of the benchmark
data(trichoptera); data(barents); data(oaks); data(mollusk); data(microcosm)
tri <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)
mol <- suppressMessages(prepare_data(mollusk$Abundance, mollusk$Covariate))
case <- function(name, data, formula, zi) list(name = name, data = data, formula = formula, zi = zi)
cases <- list(
  case("microcosm site col", microcosm, Abundance ~ 0 + site + offset(log(Offset)), "col"),
  case("microcosm site|site", microcosm, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"),
  case("trichoptera single", tri, Abundance ~ 1 + offset(log(Offset)), "single"),
  case("trichoptera col", tri, Abundance ~ 1 + offset(log(Offset)), "col"),
  case("barents single", barents, Abundance ~ 1 + offset(log(Offset)), "single"),
  case("barents col", barents, Abundance ~ 1 + offset(log(Offset)), "col"),
  case("mollusk single", mol, Abundance ~ 1 + offset(log(Offset)), "single"),
  case("mollusk col", mol, Abundance ~ 1 + offset(log(Offset)), "col"),
  case("mollusk site|site", mol, Abundance ~ 0 + site + offset(log(Offset)) | 0 + site, "covar"),
  case("oaks single", oaks, Abundance ~ 1 + offset(log(Offset)), "single"),
  case("oaks col", oaks, Abundance ~ 1 + offset(log(Offset)), "col"),
  case("oaks tree col", oaks, Abundance ~ 0 + tree + offset(log(Offset)), "col"),
  case("oaks tree|tree", oaks, Abundance ~ 0 + tree + offset(log(Offset)) | 0 + tree, "covar"))
