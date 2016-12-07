Data Structure

Package for parsing netMHC predictions. Data is stored in a


General Info:
 
 
    class PredictionCollection
 
         PredictionCollection.dictionary_collections
                
                [ { protein_1 : {'High Affinity Regions : [region1, region 2, ...],
                                 'ProtSeq : 'SNFEDNSKDAJSNFWIAIAIJFNWKSAM....',
                                 'Predictions : [ Prediction,
                                                  Prediction,
                                                  ...
                                                  Prediction.Swaps.swap
                                                
                                                    [ { swap : Prediction },
                                                      { swap : Prediction },
                                                        ...,
                                                      { swap : Prediction }
                                                    ]
                                                 ],       
                   
                  { protein_3 : {'High Affinity Regions : [region1, region 2, ...],
                                 'ProtSeq : 'SNFEDNSKDAJSNFWIAIAIJFNWKSAM....',
                                 'Predictions : [ Prediction,
                                                  Prediction,
                                                  ...
                                                  Prediction.Swaps.swap
                                                
                                                    [ { swap : Prediction },
                                                      { swap : Prediction },
                                                        ...,
                                                      { swap : Prediction }
                                                    ]
                                                 ],     
                                 
                    ...,                   
                                    
                  { protein_1 : {'High Affinity Regions : [region1, region 2, ...],
                                 'ProtSeq : 'SNFEDNSKDAJSNFWIAIAIJFNWKSAM....',
                                 'Predictions : [ Prediction,
                                                  Prediction,
                                                  ...
                                                  Prediction.Swaps.swap
                                                
                                                    [ { swap : Prediction },
                                                      { swap : Prediction },
                                                        ...,
                                                      { swap : Prediction }
                                                    ]
                                                 ],     `
         
                             
More specific information about the Prediction class:                                                                    


    class Prediction #Holds all the info about a single prediction.
    
        #I.e.:
        Prediction.allele
        Prediction.original_position
        Prediction.protein
        Prediction.affinity_level
        Prediction.rank
        Prediction.core_pep
        Prediction.h_avg_ranks
        Prediction.n_binders 
        Prediction.allele 
        Prediction.nmer`

And as shown above, it also holds data about the possible swaps. For high-affinity peptides (user-defined thershold), the Swap class has as attribute a dictionary with the following structure:

    Prediction.Swap.swaps

        {'swap1': Prediction,
         'swap2': Prediction,
         'swap3': Prediction,
          ...
         'swapN': Prediction}

     
Where `swapK` is a singly-modified peptide from the original, high-affinity peptide. 
        

                    
                                        
