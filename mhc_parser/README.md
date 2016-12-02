Data Structure

Package for parsing netMHC predictions. Data is stored in a



   ```
    class PredictionCollection
    
        PredictionCollection.dictionary_collections
            
            [ { protein_1 : [ Prediction
                              Prediction
                              ...
                              Prediction.Swaps.swap
                            
                                [ { swap : Prediction },
                                  { swap : Prediction },
                                    ...,
                                  { swap : Prediction }
                                ]
                             ],       
               
              { protein_2 : [ Prediction
                              Prediction
                              ...,
                              Prediction.Swaps.swap
                        
                                [ { swap : Prediction },
                                  { swap : Prediction },
                                  ...,
                                  { swap : Prediction }
                                ]
                                  
                             ],
                             
                ...,                   
                                
              { protein_N : [ Prediction
                              Prediction
                              ...,
                              Prediction.Swaps.swap
                        
                                [ { swap : Prediction },
                                  { swap : Prediction },
                                    ...
                                  { swap : Prediction }  
                                ]
                                  
                             ]
                            
            ] 
                             
                                                                         
                            
                              
                                        
