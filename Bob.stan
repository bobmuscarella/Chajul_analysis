data {
    int<lower=0> nspecies;
    int<lower=0> nindiv;
    int<lower=0> nplot;
    int<lower=0> n;
    int<lower=0> species[n];
    int<lower=0> indiv[n];
    int<lower=0> plot[n];
    vector[n] growth;
    vector[n] plotagb;
    vector[n] dbh;
    vector[nspecies] trait;
}


parameters {
    real beta_1[nspecies];
    real beta_2[nspecies];
    real beta_3[nspecies];
    real mu_beta[3];
    real beta_t[2];
    real <lower=0.00001> tau[5];
    real plot_effect[nplot];
    real indiv_effect[nindiv];


}

transformed parameters {
   real mu[n];
for(i in 1:n)
{
mu[i]<-exp(beta_1[species[i]]+beta_2[species[i]]+beta_3[species[i]]+plot_effect[plot[i]]+indiv_effect[indiv[i]]);

#mu[i]<-beta_1[species[i]]+beta_2[species[i]]+beta_3[species[i]]+plot_effect[plot[i]]+indiv_effect[indiv[i]]; # or without exp link
}
}

model {

for(i in 1:5)
{
tau[i]~gamma(100,100);
}

for(i in 1:3)
{
mu_beta[i]~normal(0,100);
}

for(i in 1:2)
{
beta_t[i]~normal(0,100);
}

    for(i in 1:nplot)  
{
plot_effect[i]~normal(0,tau[5]);
}

    for(i in 1:nindiv)  
{
indiv_effect[i]~normal(0,tau[4]);
}
    
  for(i in 1: nspecies)
    {
        beta_1[i]~normal(mu_beta[1]+beta_t[1]*trait[i],tau[2]);  
        beta_2[i]~normal(mu_beta[2]+beta_t[2]*trait[i],tau[3]);  
        beta_3[i]~normal(mu_beta[3],tau[4]);  
    }

    for(i in 1:n)  
    {   
      growth[i]~normal(mu[i],tau[1]);  
    }
}

generated quantities { 
vector[n] ynew;
for(i in 1:n)
    {
ynew[i]<-normal_rng(mu[i],tau[1]);  # then you can sample ynew and it’s the fitted y
}
}


