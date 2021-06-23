#####################################################################################################
# Function for FDR corrected hypergeometric enrichment test
#####################################################################################################

# snow and rlecuyer packages should be installed


## Arguments of HyperGeomFDR function:
#	steps:		the rounds of simulations (a single number)
#	pool:		background genes (character vector)
#	select:		genes to investigate (character vector)
#	DB: 		the genes set used for enrichment analysis (character list)
#	nthreads:	number of threads to use (a single number)

## Description of the hypergeometric test
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# Arguments:
#        q: vector of quantiles representing the number of white balls
#           drawn without replacement from an urn which contains both
#           black and white balls.
#        m: the number of white balls in the urn.
#        n: the number of black balls in the urn.
#        k: the number of balls drawn from the urn.
# 
# x=length(intersect(select,DB_i))    	#Number of common genes between DB and select
# m=length(intersect(pool,DB_i))        #Number of common genes between DB and pool
# n=length(pool)-length(intersect(pool,DB_i))     #Number of non-pool genes among DB (setdiff)
# k=length(select)                    	#Number of genes in select
# P_val=dhyper(length(intersect(select,DB_i)), length(intersect(pool,DB_i)), length(pool)-length(intersect(pool,DB_i)), length(select))
# 
# wikipedia
# N = length(pool)
# K = length(intersect(pool,DB_i))
# n = length(select)
# k = length(intersect(select,DB_i))
# 
# (choose(length(intersect(pool,DB_i)),length(intersect(select,DB_i)))*
#      choose(length(pool)-length(intersect(pool,DB_i)),length(select)-length(intersect(select,DB_i))))/
#      choose(length(pool),length(select))

# m= intersect(select,DB_i)
# n= intersect(bg,DB_i)-m
# k=length (DB_i)


# P_val=(choose(length(intersect(BG,DB_i)),length(intersect(IG,DB_i)))*choose(length(BG)-length(intersect(BG,DB_i)),length(IG)-length(intersect(IG,DB_i))))/choose(length(BG),length(IG))

HyperGeomFDR=function(steps, pool, select, DB, nthreads=4) {
#	steps=3000
#	pool=row.names(fc.mx)
#	select=my.selection
#	DB=DB_genes
#	nthreads=4
#save(steps,pool,select,DB,nthreads, file="tmp1.Rdata")
#load(file="tmp1.Rdata")	
	
    DB_names=names(DB)
    number_of_categories_in_DB=length(DB)
    size_of_pool=length(pool)
    size_of_select=length(select)
    
    DB_in_select=integer(number_of_categories_in_DB)  #I pre-allocate arrays for the results #  elore lefoglalok ures tomboket az eremények számára
    DB_in_pool=integer(number_of_categories_in_DB)
    Genes_in_DB=integer(number_of_categories_in_DB)
    P_val=double(number_of_categories_in_DB)
    R_obs=integer(number_of_categories_in_DB)
    
    # for every DB entity in the DB list
    for (i in 1:number_of_categories_in_DB) {
        # create a vector of genes connected to the i-th DB category
        DB_i=DB[[i]]
        # hypergometric test
        DB_in_select[i]=length(intersect(select,DB_i))	#q: number of common genes between a DBterm and select
        DB_in_pool[i]=length(intersect(pool,DB_i))	#m: number of common genes between DBterm and BackGround
        Genes_in_DB[i]=length(DB_i)			
							#n:  number of non-pool genes among DB
							#k: number of genes in select
	# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        P_val[i]=1-phyper(DB_in_select[i]-1, DB_in_pool[i], size_of_pool-DB_in_pool[i], size_of_select)  ## TODO itt a size_select kell?
    }
	
	P_val_round=round(P_val, digits=15) ## can change the digits, this is important for the precision of '0' is R
    for (i in 1:number_of_categories_in_DB) { # TODO ez egy rang szamitas.
    	 R_obs[i]=sum(P_val_round<=P_val_round[i])
    }
    P_val_df=data.frame(DB_names, DB_in_select, DB_in_pool, Genes_in_DB, P=P_val, P_adj_Bonf=p.adjust(P_val, method="bonferroni"), P_adj_BH=p.adjust(P_val, method="BH"), R_obs)
    
    ######
    # simualtion
    ######
    R_exp=integer(number_of_categories_in_DB)
    # random sampling from pool (background genes)
    # The time consuming step. The simulation here can be parallelized
   	require(snow)
	require(rlecuyer)
	
	seeds=sample(seq(1e4,1e6),6) # max number of seeds for RNGstream is 6
	cl=makeCluster(nthreads, type="SOCK")
	clusterSetupRNG(cl, type='RNGstream', seed=seeds)
	P_Sim_vec=clusterApply(cl, rep(ceiling(steps/nthreads), nthreads), sim_hyperGeom, pool, select,DB, P_val_df$DB_in_pool ) # return a list
	stopCluster(cl)
	
	
	
	P_Sim_vec=as.vector(unlist(P_Sim_vec))
	P_Sim_round=round(P_Sim_vec, digits=15)
	P_Sim_round=sort(P_Sim_round)
	cnt.of.ones<-sum(P_Sim_round==1)
	for (l in 1:length(P_val_df$P)) { # TODO ezt is betenni a parhuzamos szamitasba
	    if(P_val_round[l]>=1)
		{
			R_exp[l]= length(P_Sim_round) # ez a resz gyorsit sokat
		}else{
		 # R_exp[l]=sum(P_Sim_round<=P_val_round[l])
		 # ezt a kifejezest lecserelem egy remelhetoleg gyorsabb megoldasra: binaris keresesre
			
			t<-P_val_round[l]
			
			a1<-1
			a2<-length(P_Sim_round)-cnt.of.ones
			#P_Sim_round[a2+1] # ez mar 1
			while(a1 !=a2 & a1+1 !=a2) # TODO egyenlo lehet egyaltalan?
			{	
				a3<-ceiling((a1+a2)/2)
				if( P_Sim_round[a3]<=t)
				{
					a1<-a3
				}else
				{
					a2<-a3	
				}
		 		R_exp[l]<-a2
	    	}
		}
	}
	P_val_df$R_exp=R_exp/steps
    P_val_df$FDR=P_val_df$R_exp/R_obs
    return(P_val_df)
}

sim_hyperGeom=function(steps, pool, select, DB,DB_in_pool) {
    DB_names=names(DB)
    num_DB=length(DB_names)
    P_Sim_mat=matrix(numeric(num_DB*steps), ncol=steps)
    size_pool=length(pool)
    size_select=length(select)

    for (j in 1:steps) {
        Rand.select=sample(pool, size_select)       
        for (i in 1:num_DB) {
            DB_i=DB[[i]]
            # hypergometric test
            P_Sim_mat[i,j]=1-phyper(length(intersect(Rand.select, DB_i))-1, DB_in_pool[i], size_pool-DB_in_pool[i], size_select)
        }
    }
    return(as.vector(P_Sim_mat))
}


