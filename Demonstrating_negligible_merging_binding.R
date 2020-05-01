#### Carla Bardua ####
#### This code is provided as an example of how to create negligible regions, and then 
#### how to combine this negligible region data with your present region data.
#### finally, I demonstrate how to combine all individual bone patches into one dataset
#### I am sure there will be a neater way to do this, but I am sharing this in the hope that 
#### it may be of some help!

#### This code assumes you already have your raw landmark and curve data saved as a 3D array
#### here, my landmark and curve data is called "subsampled.lm" and has these dimensions:
dim(subsampled.lm)
[1] 640   3 174
#### So in this example, I have 640 landmarks/curves, the data are 3D, and I have 174 specimens

##### creating negligible regions ####
## first you need a spreadsheet detailing which specimens have which negligible regions
Region_info=read.csv(file="C:/Users/carlb/Box Sync/Carla Bardua/Frogs/Region_information.csv", header=TRUE) 
## details of which specimens have absent region- put 'NA' when region is absent for each specimen

SphV=Region_info$SPH.V ## look at the relevant column- here I am looking at the SphV (Sphenethmoid ventral) region
to_remove=na.omit(SphV) ## look at which rows have NA values- so which specimens dont have a SphV region
to_remove ## see here for list of specimens that are NA- here it was specimens 57,62,79,143
remove_SphV=c(57,62,79,143) ## create vector of these specimens

### now you take your landmark/curve data (as a 3D array), which I have called "subsampled.lm', and you divide 
### this into 2 datasets- one with all specimens that have the SphV region, and one with all 
### specimens that do NOT have a sphV region.
subsampled.lm_SphV_present=subsampled.lm[,,-remove_SphV] ## this you will then patch as normal
subsampled.lm_SphV_absent=subsampled.lm[,,remove_SphV] ## these specimens need the negligible region

#################################
#### choose a landmark position that best represents where this negligible region should go
#### e.g. my SphV region is closest to the parasphenoid, so I chose a landmark on the parasphenoid
substitute_LM10=subsampled.lm_SphV_absent[10,,] ##  want negligible region on Parasphenoid tip which is LM10
dim(substitute_LM10) ## 4 specs

N=41 # I want 41 patch landmarks- since this was the number of patch points this region will have when present

### now, for each specimen, take its "landmark 10" and replicate this coordinate 41 times
SphV_absent_patch=array(dim=c(N, 3, dim(subsampled.lm_SphV_absent)[[3]]))
dimnames(SphV_absent_patch)[3]=dimnames(subsampled.lm_SphV_absent)[3]

i=1
for (i in 1:N)
{
  SphV_absent_patch[i,,] <- substitute_LM10
  i=i+1
}

dim(SphV_absent_patch)
SphV_absent_patch_N_4=SphV_absent_patch
save(SphV_absent_patch_N_4, file="./patched_datasets/SphV_absent_patch_N_4.R")

## All of my patches were saved as separate files like this, so at the end I load them all into R and combine
## them into one array

##### COMBINING PRESENT AND NEGLIGIBLE DATA FROM EACH REGION
#### if you've patched a region in batches (e.g. you may want to patch 50 specimens at a time for ease),
#### combine them all here too
## So for the SphV region I patched the "present" specimens in 2 batches, and the "absent" specimens in 2 batches for this bone

#### SphV ####
#### load the patched datasets (that you may have patched in batches)
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Final/SphV_patch_ordered_Fl_KD_NHM_N_108.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Final/SphV_absent_patch_N_46.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Final/Ornam_etc_SphV_patch.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Final/Ornam_SphV_absent_patch_N_7.R")

patched_collated=abind::abind(SphV_patch_ordered_Fl_KD_NHM_N_108,
                              SphV_absent_patch_N_46,
                              Ornam_etc_SphV_patch,
                              Ornam_SphV_absent_patch_N_7,
                              along = 3)
## here you combine all specimens for this bone, including both present and negligible regions

## make in alphabetical order as order here will be in batch order
N = dim(patched_collated)[[1]]
data_array <- two.d.array(patched_collated)
data_array_ordered= data_array[order(rownames(data_array)),]
patched_collated_ordered <- arrayspecs(data_array_ordered, N, 3)
dimnames(patched_collated_ordered)[[3]]
spheres3d(patched_collated_ordered[,,2], radius = 0.05)
dim(patched_collated_ordered)

SphV_patch_final=patched_collated_ordered

#save(SphV_patch_final, file="E:/LOOK HERE/Frogs/PATCHING/patched_datasets/SphV_patch_final.R")

####### COMBINING DIFFERENT REGIONS INTO ONE ARRAY ####
#### once you've created your final patched datasets for each bone (combining present and negligible regions),
#### you can then load in each patched bone and combine these into one dataset

## they should all be in alphabetical order now
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Merged_final/Premax_dorsal_patch_final.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Merged_final/Premax_ventral_patch_final.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Merged_final/PS_patch_final.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Merged_final/Pt_patch_final_2.R")
load("E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Merged_final/MaxV_patch_final.R")
## etc

Final_patch_combined=abind::abind(Premax_dorsal_patch_final,Premax_ventral_patch_final,
                                    PS_patch_final,Pt_patch_final_2,MaxV_patch_final,
                                    along=1)

dim(Final_patch_combined)
dimnames(Final_patch_combined)

save(Final_patch_combined, file="E:/LOOK HERE/Frogs/PATCHING/patched_datasets/Final_patch_combined.R")

## this is your final patch data which you can then bind with your landmarks and curve data