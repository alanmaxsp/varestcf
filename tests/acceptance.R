# One deterministic acceptance suite; no BLUPF90, private files or compiled code.
library(varestcf)
source(system.file("examples", "synthetic_models.R", package = "varestcf"))
checks <- 0L
assert <- function(value) {
  if (!isTRUE(value)) stop("Acceptance assertion failed", call. = FALSE)
  checks <<- checks + 1L
}
close <- function(a, b, tol = 1e-9) assert(max(abs(a - b)) <= tol * max(1, abs(b)))
reject <- function(expr, text) {
  message <- tryCatch({force(expr); "NO ERROR"}, error = function(e) conditionMessage(e))
  assert(grepl(text, message, fixed = TRUE))
}
run <- function(model, ...) do.call(reconstruct_focal_pev,
                                  c(model[setdiff(names(model), "y")], list(...)))
models <- synthetic_models();results <- list();oracles <- list()
for (name in names(models)) {
  model <- models[[name]];oracle <- dense_reference(model)
  got <- run(model, block_size = 3)
  close(got$Sigma, oracle$Sigma)
  assert(identical(rownames(got$Sigma), model$focal_ids))
  assert(max(got$diagnostics$solves$backward_error) < 1e-10)
  assert(all(diag(got$Sigma) > 0))
  # Independent explicit centering-matrix formula, including fixed-effect uncertainty.
  n <- length(oracle$ebv);P <- diag(n) - matrix(1 / n, n, n)
  expected <- as.numeric(crossprod(oracle$ebv, P %*% oracle$ebv) + sum(diag(P %*% oracle$Sigma))) / n
  cf <- varest_cf(oracle$ebv, got$Sigma)
  close(cf$estimate, expected)
  close(varest_cf(oracle$ebv, got$Sigma, "n-1")$estimate, expected * n / (n - 1))
  close(varest_cf(oracle$ebv + 100, got$Sigma)$estimate, cf$estimate)
  order <- rev(seq_len(n))
  close(varest_cf(oracle$ebv[order], got$Sigma)$estimate, cf$estimate)
  results[[name]] <- got;oracles[[name]] <- oracle
  cat("PASS model:", name, "\n")
}

# Reordering the precision independently is aligned by coefficient identifiers.
model <- models$animal_combined_relationship;perm <- rev(seq_len(ncol(model$Z)))
model$random_precision <- model$random_precision[perm, perm]
close(run(model, block_size = 1)$Sigma, results$animal_combined_relationship$Sigma)
model$Z <- model$Z[, perm]
close(run(model, block_size = 20)$Sigma, results$animal_combined_relationship$Sigma)
# General correlated residual precision and its named record permutation.
model <- models$animal_pedigree_relationship;n <- nrow(model$Z)
W <- diag(.2, n) + matrix(.002, n, n)
dimnames(W) <- list(rownames(model$Z), rownames(model$Z))
model$residual_precision <- W
oracle <- dense_reference(model);got <- run(model)
close(got$Sigma, oracle$Sigma)
model$residual_precision <- W[n:1, n:1]
close(run(model)$Sigma, got$Sigma)
# Models without fixed effects are supported; duplicate fixed columns are rejected.
model <- models$animal_pedigree_relationship;model$X <- model$X[, FALSE, drop = FALSE]
close(run(model)$Sigma, dense_reference(model)$Sigma)
model <- models$animal_pedigree_relationship;model$X <- cbind(model$X, model$X[, 1])
reject(run(model), "full column rank")

# A Schur complement independently checks the uncertainty of the fixed effects.
o <- oracles$animal_combined_relationship;p <- o$p;r <- nrow(o$C) - p
B <- o$C[seq_len(p), seq_len(p)];F <- o$C[p + seq_len(r), seq_len(p)]
Di <- solve(o$C[p + seq_len(r), p + seq_len(r)])
schur <- Di + Di %*% F %*% solve(B - t(F) %*% Di %*% F) %*% t(F) %*% Di
idx <- match(models$animal_combined_relationship$focal_ids, colnames(models$animal_combined_relationship$Z))
close(schur[idx, idx], o$Sigma)
assert(max(abs(Di[idx, idx] - o$Sigma)) > 1e-4)

# random regression projection must retain between-animal cross-coefficient covariances.
for (name in c("random_regression_homogeneous", "random_regression_heterogeneous")) {
  o <- oracles[[name]];S <- results[[name]]$Sigma;n <- nrow(S) / 2
  L <- cbind(diag(n), .7 * diag(n));projected <- L %*% S %*% t(L)
  labels <- paste0("focal", seq_len(n));dimnames(projected) <- list(labels, labels)
  u <- setNames(as.numeric(L %*% o$ebv), labels)
  close(projected, L %*% o$Sigma %*% t(L))
  bad <- S;cross <- S[1:n, n + 1:n];bad[1:n, n + 1:n] <- diag(diag(cross))
  bad[n + 1:n, 1:n] <- diag(diag(cross))
  bad <- L %*% bad %*% t(L);dimnames(bad) <- list(labels, labels)
  close(diag(bad), diag(projected))
  assert(abs(varest_cf(u, projected)$estimate - varest_cf(u, bad)$estimate) > 1e-5)
}

# Input contracts: names, covariance completeness, scale and resource guards.
o <- oracles$animal_pedigree_relationship;S <- o$Sigma;u <- o$ebv
reject(varest_cf(unname(u), S), "ebv names")
bad <- S;bad[1, 2] <- NA_real_;reject(varest_cf(u, bad), "nonfinite")
bad <- S;bad[1, 2] <- bad[1, 2] + .2;reject(varest_cf(u, bad), "symmetric")
bad <- S;bad[1, 1] <- -1;reject(varest_cf(u, bad), "semidefinite")
bad <- S;rownames(bad)[1] <- "unknown";reject(varest_cf(u, bad), "identifiers")
reject(varest_cf(u, Matrix::Matrix(S, sparse = TRUE)), "dense matrix")
zero <- S * 0;close(varest_cf(u, zero)$estimate, sum((u - mean(u))^2) / length(u))
common <- S * 0 + 1;close(varest_cf(u, common)$estimate, varest_cf(u, zero)$estimate)
reject(varest_cf(u[1], S[1, 1, drop = FALSE]), "at least two")
model <- models$animal_pedigree_relationship;model$focal_ids <- rep(model$focal_ids[1], 2)
reject(run(model), "unique")
model <- models$animal_pedigree_relationship;model$focal_ids[1] <- "unknown";reject(run(model), "present in Z")
model <- models$animal_pedigree_relationship;model$residual_precision <- -1;reject(run(model), "positive")
model <- models$animal_pedigree_relationship;model$random_precision[1, 1] <- -1;reject(run(model), "positive definite")
model <- models$animal_pedigree_relationship;rownames(model$X) <- rev(rownames(model$X));reject(run(model), "record orders")
reject(run(models$animal_pedigree_relationship, block_size = 0), "block_size")
reject(run(models$animal_pedigree_relationship, memory_budget_gib = 1e-12), "memory_budget_gib")
reject(run(models$animal_pedigree_relationship, tolerance = NA_real_), "tolerance")
# High-level formula workflow: compared with the independent matrix fixtures above.
ped <- data.frame(animal=paste0("animal",1:8),
 sire=c(0,0,1,1,3,3,5,5), dam=c(0,0,2,2,4,4,6,6))
for(j in 2:3) ped[[j]] <- ifelse(ped[[j]]==0, NA_character_, paste0("animal",ped[[j]]))
dat <- data.frame(animal=rep(ped$animal,each=3), group=rep(c(-.5,.5),12),
 time=seq(-1.2,1.2,length.out=24), weight=as.numeric(models$animal_pedigree_relationship$y))
base_args <- list(data=dat, formula=weight~group, animal="animal", relationship_inverse=models$animal_pedigree_relationship$random_precision*4,
 focal_ids=models$animal_pedigree_relationship$focal_ids, genetic_variance=4, residual_variance=5)
high <- do.call(reconstruct_focal_pev,base_args)
close(high$Sigma,oracles$animal_pedigree_relationship$Sigma);close(high$ebv,oracles$animal_pedigree_relationship$ebv)
close(varest_cf(high)$estimate,varest_cf(oracles$animal_pedigree_relationship$ebv,oracles$animal_pedigree_relationship$Sigma)$estimate)
assert(inherits(high,"varestcf_fit"));assert(inherits(varest_cf(high),"varestcf_result"))
# A supplied H or its inverse avoids pedigree and manual joint-precision assembly.
arg <- base_args;arg$relationship_inverse <- NULL;arg$relationship <- solve(models$animal_combined_relationship$random_precision*4)
hs <- do.call(reconstruct_focal_pev,arg)
close(hs$Sigma,oracles$animal_combined_relationship$Sigma);close(hs$ebv,oracles$animal_combined_relationship$ebv)
arg$relationship_inverse <- solve(arg$relationship);arg$relationship <- NULL
close(do.call(reconstruct_focal_pev,arg)$Sigma,hs$Sigma)
# Redundant predictors and a constant factor are handled by identifying their column space.
arg <- base_args;arg$data$duplicate <- dat$group;arg$data$constant <- factor(rep("one",24))
arg$formula <- weight~group+duplicate+constant
close(do.call(reconstruct_focal_pev,arg)$Sigma,high$Sigma)
# Heterogeneous random regression: projections and EBV are prepared without caller-built matrices.
arg <- base_args;arg$relationship_inverse <- models$animal_combined_relationship$random_precision*4
arg$focal_ids <- paste0("animal",c(3,5,6,8));arg$random <- ~1+time
arg$genetic_variance <- matrix(c(4,.3,.3,.7),2)
arg$residual_variance <- ~exp(log(5)+.25*time)
hr <- do.call(reconstruct_focal_pev,arg)
close(hr$Sigma,oracles$random_regression_heterogeneous$Sigma);close(hr$ebv,oracles$random_regression_heterogeneous$ebv)
grid <- c(-.7,0,.7);automatic <- varest_cf(hr,at=grid)
for(i in seq_along(grid)){
 n <- 4;L <- cbind(diag(n),grid[i]*diag(n));S <- L%*%oracles$random_regression_heterogeneous$Sigma%*%t(L)
 dimnames(S) <- list(arg$focal_ids,arg$focal_ids)
 u <- setNames(as.numeric(L%*%oracles$random_regression_heterogeneous$ebv),arg$focal_ids)
 close(automatic$estimate[i],varest_cf(u,S)$estimate)
}
reject(varest_cf(hr),"at=")
close(varest_cf(hr,at=data.frame(time=grid))$estimate,automatic$estimate)
# Multivariate response, including missing traits and correlated residuals.
arg <- base_args;arg$relationship_inverse <- models$animal_combined_relationship$random_precision*4
arg$focal_ids <- paste0("animal",c(3,5,6,8));arg$formula <- cbind(t1,t2,t3,t4)~1
for(j in 1:4) arg$data[[paste0("t",j)]] <- ifelse(rep(1:4,6)==j,dat$weight,NA_real_)
arg$genetic_variance <- diag(c(4,5,3,6))+.3;arg$residual_variance <- c(5,6,7,8)
multi <- do.call(reconstruct_focal_pev,arg)
close(multi$Sigma,oracles$multivariate_animal$Sigma);close(multi$ebv,oracles$multivariate_animal$ebv)
assert(nrow(varest_cf(multi)$table)==4L)
arg$formula <- cbind(t1,t2)~group;arg$data$t1 <- dat$weight;arg$data$t2 <- dat$weight+cos(1:24)
arg$data$t2[c(2,7,10)] <- NA_real_;arg$genetic_variance <- matrix(c(4,.3,.3,2),2)
arg$residual_variance <- matrix(c(5,.4,.4,6),2)
two <- do.call(reconstruct_focal_pev,arg)
# Independent stack and residual marginalization, not a subblock of a full inverse.
present <- !is.na(as.matrix(arg$data[,c("t1","t2")]));loc <- which(present,arr.ind=TRUE)
XX <- matrix(0,nrow(loc),4);ZZ <- matrix(0,nrow(loc),16);WW <- matrix(0,nrow(loc),nrow(loc))
for(i in seq_len(nrow(loc))){
 rr <- loc[i,1];tt <- loc[i,2];XX[i,(tt-1)*2+1:2] <- c(1,dat$group[rr])
 ZZ[i,(tt-1)*8+match(dat$animal[rr],ped$animal)] <- 1
}
for(rr in unique(loc[,1])){ii <- which(loc[,1]==rr);WW[ii,ii] <- solve(arg$residual_variance[loc[ii,2],loc[ii,2],drop=FALSE])}
DD <- cbind(XX,ZZ);CC <- crossprod(DD,WW%*%DD)
CC[4+1:16,4+1:16] <- CC[4+1:16,4+1:16]+kronecker(solve(arg$genetic_variance),arg$relationship_inverse)
ii <- 4+c(c(3,5,6,8),8+c(3,5,6,8));close(two$Sigma,solve(CC)[ii,ii])
# User-facing errors and response handling.
arg <- base_args;arg$data$weight[1] <- NA_real_;partial <- do.call(reconstruct_focal_pev,arg)
assert(partial$model$excluded_response_values==1L)
arg <- base_args;arg$data$group[1] <- NA_real_;reject(do.call(reconstruct_focal_pev,arg),"missing values")
arg <- base_args;arg$relationship <- diag(8);reject(do.call(reconstruct_focal_pev,arg),"exactly one")
arg <- base_args;arg$genetic_variance <- NULL;reject(do.call(reconstruct_focal_pev,arg),"Supply genetic_variance")
arg <- base_args;arg$focal_ids <- "not-an-animal";reject(do.call(reconstruct_focal_pev,arg),"absent")
arg <- base_args;arg$residual_variance <- -1;reject(do.call(reconstruct_focal_pev,arg),"positive")
reject(varest_cf(hr,at=grid,component="(Intercept)"),"Choose at or component")
# A class-specific residual variance is equivalent to its explicit row vector.
arg <- base_args;arg$random <- ~1+time
arg$genetic_variance <- matrix(c(4,.3,.3,.7),2)
arg$data$residual_class <- rep(1:3,each=8)
arg$residual_variance <- c(2,5,9)[arg$data$residual_class]
class_fit <- do.call(reconstruct_focal_pev,arg)
arg$residual_variance <- ~c(2,5,9)[residual_class]
close(do.call(reconstruct_focal_pev,arg)$Sigma,class_fit$Sigma)
arg$residual_variance <- 5
homogeneous_fit <- do.call(reconstruct_focal_pev,arg)
assert(max(abs(class_fit$Sigma-homogeneous_fit$Sigma))>1e-5)
assert(max(abs(class_fit$ebv-homogeneous_fit$ebv))>1e-5)
# Nuisance terms affect genetic uncertainty and must remain in the equations.
assert(max(abs(results$repeatability$Sigma-results$animal_pedigree_relationship$Sigma))>1e-5)
for(name in c("repeatability","direct_maternal")){
 full_fit <- do.call(reconstruct_focal_pev,models[[name]])
 close(full_fit$ebv,oracles[[name]]$ebv)
}
# Automatic multiple-term construction; compare with independent dense equations.
Ai <- models$animal_pedigree_relationship$random_precision*4
Hi <- models$animal_combined_relationship$random_precision*4
focals <- models$animal_pedigree_relationship$focal_ids
effects <- list(additive=list(id="animal",relationship_inverse=Ai,variance=4),
                permanent=list(id="animal",variance=2))
repeated_args <- list(data=dat,formula=weight~group,random_effects=effects,
                     focal_effect="additive",focal_ids=focals,residual_variance=5)
repeated <- do.call(reconstruct_focal_pev,repeated_args)
close(repeated$Sigma,oracles$repeatability$Sigma)
close(repeated$ebv,oracles$repeatability$ebv)
close(varest_cf(repeated)$estimate,varest_cf(oracles$repeatability$ebv,oracles$repeatability$Sigma)$estimate)
arg <- repeated_args;arg$random_effects <- rev(effects)
close(do.call(reconstruct_focal_pev,arg)$Sigma,repeated$Sigma)
arg$random_effects$permanent$levels <- rev(rownames(Ai))
close(do.call(reconstruct_focal_pev,arg)$Sigma,repeated$Sigma)
arg <- repeated_args;arg$random_effects <- effects["additive"]
close(do.call(reconstruct_focal_pev,arg)$Sigma,oracles$animal_pedigree_relationship$Sigma)
# Direct/maternal covariance: the missing IDs explicitly mean zero incidence.
maternal_data <- dat;maternal_data$dam <- rep(ped$dam,each=3)
dm <- list(direct=list(id="animal",relationship_inverse=Ai),
           maternal=list(id="dam",relationship_inverse=Ai,allow_missing=TRUE))
Gdm <- matrix(c(4,.4,.4,2),2,dimnames=list(c("direct","maternal"),c("direct","maternal")))
maternal_args <- list(data=maternal_data,formula=weight~group,random_effects=dm,
 focal_effect="direct",focal_ids=focals,residual_variance=5,
 correlated_effects=list(list(effects=c("direct","maternal"),covariance=Gdm)))
direct_fit <- do.call(reconstruct_focal_pev,maternal_args)
close(direct_fit$Sigma,oracles$direct_maternal$Sigma);close(direct_fit$ebv,oracles$direct_maternal$ebv)
arg <- maternal_args;arg$random_effects$maternal$relationship_inverse <- Ai[8:1,8:1]
arg$correlated_effects[[1]]$effects <- c("maternal","direct")
close(do.call(reconstruct_focal_pev,arg)$Sigma,direct_fit$Sigma)
arg <- maternal_args;arg$focal_effect <- "maternal"
maternal_oracle <- models$direct_maternal;maternal_oracle$focal_ids <- paste0("maternal:",focals)
close(do.call(reconstruct_focal_pev,arg)$Sigma,dense_reference(maternal_oracle)$Sigma)
assert(direct_fit$model$random_effects$maternal$zero_incidence_rows==6L)
# Independent genetic/permanent random regressions with heterogeneous residuals.
rr <- models$random_regression_heterogeneous
PE <- cbind(models$animal_pedigree_relationship$Z,models$animal_pedigree_relationship$Z*dat$time)
colnames(PE) <- paste0("pe",seq_len(ncol(PE)))
Gpe <- matrix(c(2,.1,.1,.5),2)
rr_full <- rr;rr_full$Z <- cbind(rr$Z,PE)
rr_full$random_precision <- as.matrix(Matrix::bdiag(rr$random_precision,kronecker(solve(Gpe),diag(8))))
dimnames(rr_full$random_precision) <- list(colnames(rr_full$Z),colnames(rr_full$Z))
rr_args <- list(data=dat,formula=weight~group,
 random_effects=list(additive=list(id="animal",random=~1+time,relationship_inverse=Hi,
                                   variance=matrix(c(4,.3,.3,.7),2)),
                     permanent=list(id="animal",random=~1+time,variance=Gpe)),
 focal_effect="additive",focal_ids=paste0("animal",c(3,5,6,8)),
 residual_variance=~exp(log(5)+.25*time))
rr_fit <- do.call(reconstruct_focal_pev,rr_args);rr_oracle <- dense_reference(rr_full)
close(rr_fit$Sigma,rr_oracle$Sigma);close(rr_fit$ebv,rr_oracle$ebv)
L <- cbind(diag(4),.7*diag(4));S <- L%*%rr_oracle$Sigma%*%t(L)
dimnames(S) <- list(rr_args$focal_ids,rr_args$focal_ids)
u <- setNames(as.numeric(L%*%rr_oracle$ebv),rr_args$focal_ids)
close(varest_cf(rr_fit,at=.7)$estimate,varest_cf(u,S)$estimate)
# Two responses and an additional permanent effect; correlated residuals, missing values.
md <- dat;md$t1 <- dat$weight;md$t2 <- dat$weight+cos(1:24);md$t2[c(2,7,10)] <- NA_real_
Gg <- matrix(c(4,.3,.3,2),2);Gp <- matrix(c(2,.1,.1,1),2);R0 <- matrix(c(5,.4,.4,6),2)
multi_args <- list(data=md,formula=cbind(t1,t2)~group,
 random_effects=list(additive=list(id="animal",relationship_inverse=Hi,variance=Gg),
                     permanent=list(id="animal",variance=Gp)),
 focal_effect="additive",focal_ids=rr_args$focal_ids,residual_variance=R0)
multi_fit <- do.call(reconstruct_focal_pev,multi_args)
ZZ_full <- cbind(ZZ,ZZ);DD_full <- cbind(XX,ZZ_full)
QQ_full <- as.matrix(Matrix::bdiag(kronecker(solve(Gg),Hi),kronecker(solve(Gp),diag(8))))
CC_full <- crossprod(DD_full,WW%*%DD_full)
CC_full[4+1:32,4+1:32] <- CC_full[4+1:32,4+1:32]+QQ_full
ii <- 4+c(c(3,5,6,8),8+c(3,5,6,8))
close(multi_fit$Sigma,solve(CC_full)[ii,ii]);assert(nrow(varest_cf(multi_fit)$table)==2L)
# Independent common-group effect with unequal numbers of levels.
gd <- dat;gd$litter <- rep(c("l1","l2","l3"),8)
arg <- repeated_args;arg$data <- gd;arg$random_effects$permanent <- list(id="litter",variance=2)
group_fit <- do.call(reconstruct_focal_pev,arg)
J <- diag(3)[match(gd$litter,c("l1","l2","l3")),,drop=FALSE]
colnames(J) <- c("l1","l2","l3");rownames(J) <- rownames(models$animal_pedigree_relationship$Z)
group_model <- models$animal_pedigree_relationship
group_model$Z <- cbind(group_model$Z,J)
group_model$random_precision <- as.matrix(Matrix::bdiag(Ai/4,diag(3)/2))
dimnames(group_model$random_precision) <- list(colnames(group_model$Z),colnames(group_model$Z))
close(group_fit$Sigma,dense_reference(group_model)$Sigma)
# Joint component covariance for two correlated random-regression terms.
rd <- dat;rd$relative <- rep(paste0("animal",c(2:8,1)),each=3)
joint_G <- diag(c(4,.7,2,.5))+.1
arg <- rr_args;arg$data <- rd
arg$random_effects <- list(direct=list(id="animal",random=~1+time,relationship_inverse=Ai),
                          related=list(id="relative",random=~1+time,relationship_inverse=Ai))
arg$focal_effect <- "direct"
arg$correlated_effects <- list(list(effects=c("direct","related"),covariance=joint_G))
joint_fit <- do.call(reconstruct_focal_pev,arg)
inc <- models$animal_pedigree_relationship$Z
inc2 <- diag(8)[match(rd$relative,rownames(Ai)),,drop=FALSE]
joint_model <- rr
joint_model$Z <- cbind(inc,inc*dat$time,inc2,inc2*dat$time)
colnames(joint_model$Z) <- paste0(rep(c("a0:","a1:","b0:","b1:"),each=8),rep(rownames(Ai),4))
rownames(joint_model$Z) <- rownames(rr$Z)
joint_model$random_precision <- kronecker(solve(joint_G),Ai)
dimnames(joint_model$random_precision) <- list(colnames(joint_model$Z),colnames(joint_model$Z))
close(joint_fit$Sigma,dense_reference(joint_model)$Sigma)
# Contracts prevent dropping terms, mixing scales and inventing cross-kernels.
arg <- repeated_args;arg$focal_effect <- NULL;reject(do.call(reconstruct_focal_pev,arg),"Select one focal_effect")
arg <- repeated_args;arg$animal <- "animal";reject(do.call(reconstruct_focal_pev,arg),"do not mix")
arg <- repeated_args;arg$random_effects$permanent$variance <- NULL
reject(do.call(reconstruct_focal_pev,arg),"Supply variance")
arg <- repeated_args;arg$random_effects$permanent$unexpected <- TRUE
reject(do.call(reconstruct_focal_pev,arg),"Invalid fields")
arg <- maternal_args;arg$random_effects$maternal$allow_missing <- FALSE
reject(do.call(reconstruct_focal_pev,arg),"allow_missing")
arg <- maternal_args;arg$random_effects$maternal$relationship_inverse <- Ai*1.1
reject(do.call(reconstruct_focal_pev,arg),"same relationship inverse")
arg <- maternal_args;arg$random_effects$direct$variance <- 99
reject(do.call(reconstruct_focal_pev,arg),"conflicts")
arg <- maternal_args;arg$correlated_effects[[1]]$covariance[1,1] <- -1
reject(do.call(reconstruct_focal_pev,arg),"positive definite")
arg <- maternal_args;arg$correlated_effects <- rep(arg$correlated_effects,2)
reject(do.call(reconstruct_focal_pev,arg),"cannot overlap")
cat("ALL ACCEPTANCE CHECKS PASSED:", checks, "\n")
