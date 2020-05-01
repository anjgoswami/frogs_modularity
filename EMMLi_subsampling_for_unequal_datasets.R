#### subsampling EMMLi ####
#### Carla Bardua ####
#### This is an example of how to run EMMLi, subsampling a dataset randomly X number of times
#### example: 
#### here I want to compare the integration patterns of frogs that have feeding vs non feeding larvae
#### my frogs with feeding larvae N = 124, whereas my frogs with non-feeding larvae N = 39
#### So I want to subsample my larger dataset (124) down to the size of the smaller dataset (39)
#### I want to do this randomly 100 times and find the average of this.

#### you will need to source this code:
source("./code from lab/newEmmliFunctions.R") 

#### here, "Y.gpa.rhs.feed" is my dataset for frogs with feeding larvae 
#### (aligned separately from the dataset of frogs without feeding larvae)
#### "modules_no_NA" is my module hypotheses for EMMLi, but no 'NA' values are allowed for this function

#### before you run this, you'll have to run the code once with one subset, in order to work out
#### what the output looks like
#### here, the output of "EMMLi_feed_random_39" has a long string of outputs inside
#### you need to extract the output that gives you the correlation values for the best supported model
#### mine looks like this:
#### EMMLi_feed_random_39$`1`$rhos_best$`By.Region.19.sep.Mod + sep.between`
#### each person's will look different though depending on the best model
#### Once you know this output, you can run the code as a loop and save this output each time

N=100 ## number of times you want to run the subsamples
Res_combined=array(dim=c(2,191,N)) ## 2 and 191 come from the dimensions of the EMMLi output

i=1
for (i in 1:N)
  
{
  EMMLi_feed_random_39 <- subSampleEMMLi(landmarks = Y.gpa.rhs.feed[,,c(sort(sample(1:124, 39, replace=FALSE)))],
                                         fractions = 1,
                                         models = modules_no_NA,
                                         min_landmark = 3, ### if too low (less than 3 I think) EMMLi won't run)
                                         return_corr = TRUE,
                                         summarise = FALSE)
  
  Res=EMMLi_feed_random_39$`1`$rhos_best$`By.Region.19.sep.Mod + sep.between`
  Res=Res[,order(colnames(Res), decreasing=FALSE)] 
  Res_combined[,,i]=Res
  i=i+1
}  

rownames(Res_combined)=rownames(Res)
colnames(Res_combined)=colnames(Res) 

Res_ave=apply(Res_combined, c(1,2), mean)

#### now plot average result, with module names and a layout if possible
plotNetwork(Res_ave, module_names = mod.names, linecolour = "darkblue", 
            title = NULL, layout = layout)

### this average result can then be compared to the result from the smaller dataset
