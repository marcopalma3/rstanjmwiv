#include pre.stan

// This Stan program defines a simple JM-WIV model. 

functions{
  #include /JMmodules/functions.stan
}

data {
  //----- Longitudinal submodels
  #include /JMmodules/data_melsm.stan
  #include /JMmodules/data_event.stan
  #include /JMmodules/data_prior.stan
}

transformed data {
  #include /JMmodules/tdata.stan
}

parameters {
  #include /JMmodules/param.stan
}

transformed parameters {
  #include /JMmodules/tparam.stan
}

model {
  #include /JMmodules/model_ALLassoc.stan
}

generated quantities {
  #include /JMmodules/gqs.stan
}
