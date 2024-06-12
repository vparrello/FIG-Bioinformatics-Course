## Projections

### 3.1 Sample Set Creation 
    * Create a sample report of hammers that lists every repgen genome inside of a sample.
    
### 3.2 Data Curation
    * What is the cut off for how many hammers are needed to classify the genome as there. Scaled to the size of the sample file
        1 MB has less hits in general than 1GB so the hammers have to account for this difference
    * Try on another sample and determine if your cut off is correct
    * Now do a whole bunch of samples. Can we use the hammer profile to characterize the sample. We do that by creating an xmatrix.
    
 ### 3.3 Creating the xmatrix   
    * Create the universal table and template for the data - helps with training data too
    * Create the x-matrix on presence and absence
    * Create a second x-matrix on scaled population

### 3.4 Machine Learning Classifier
    * Run the classifier and interpret the data
    * Run classifier on population instead of presence and absence and see which one works better.
    * Play with the parameters and see how the accuracy changes based on them
