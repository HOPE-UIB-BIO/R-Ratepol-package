fc_CharSNI = function(CharData,BandWidth) {
  
  # Code obtained and from 
  # Kelly RF, Higuera PE, Barrett CM, Feng Sheng H. 2011 A signal-to-noise index to quantify the potential for peak detection 
  #   in sediment-charcoal records. Quat. Res. 75, 11–17. (doi:10.1016/j.yqres.2010.07.011)
  
  ## Calculate signal-to-noise index (SNI) for a charcoal record.
  #
  # CharData = Matrix of input data with one row per sample, containing:
  #   Column 1: age associated with the sample (yr)
  #   Column 2: charcoal accumulation rate (CHAR) of the sample (pieces/cm^2/yr)
  #   Column 3: threshold value (pieces/cm^2/yr)
  #
  # BandWidth = Width of moving window for computing SNI 
  #
  #
  # This function computes SNI as described in Kelly et al. 2010.  Note that
  # your data must be interpolated to constant sample resolution (yr/sample)
  # before input to the function. The function makes no assumption about
  # prior analysis on the input CHAR series, i.e. any background and
  # threshold methods may be used.  However, input data should still be
  # consistent with the interpretation that a CHAR value (column 2) greater
  # than the corresponding threshold value (column 3) is a "signal" sample,
  # whereas a CHAR value below the threshold is "noise".  Refer to Kelly et
  # al. 2010 for details and discussion.
  #
  # SNI_output is a data list containing the computed SNI and related
  # data, with one row for each row in the input variable CharData:
  #   
  #   $SNI    = the SNI computed for each sample
  #   $winInd = indexes of the first and last samples included in each moving 
  #     window.  E.g. SNI_output$winInd(X) == [A, B] indicates that the 
  #     moving window used to calculate SNI for the Xth sample contained all 
  #     samples between A and B, inclusive.
  #   $popN   = the CHAR values of all samples in the noise (N) population 
  #     (samples in the moving window with CHAR below threshold)
  #   $popS   = the CHAR values of all samples in the signal (S) population 
  #     (samples in the moving window with CHAR at or above threshold)
  #   $meanN  = mean CHAR of the samples in popN
  #   $stdN   = standard deviation of the samples in popN
  #   $CF     = the "correction factor" used in computing SNI.  Equal to
  #     (v - 2)/v ,where v is the number of samples in popN
  
  # Data setup
  ages = CharData[,1];
  
  CHAR = CharData[,2];
  CHAR.mean <- mean(CHAR);
  CHAR.sd <- stats::sd(CHAR);
  CHAR <- ( CHAR - CHAR.mean ) / CHAR.sd;
  
  thresh = CharData[,3];
  thresh <- ( thresh - CHAR.mean ) / CHAR.sd;
  
  
  r = mean(diff(ages));
  
  # Preallocate some space   
  SNI_output 				= list();
  SNI_output$SNI 		= NA*numeric(length(ages));
  SNI_output$winInd	= matrix(data=NA,nrow=length(ages),ncol=2);
  SNI_output$popN 	= vector(mode="list",length=length(ages));
  SNI_output$popS 	= vector(mode="list",length=length(ages));
  SNI_output$meanN 	= NA*numeric(length(ages));
  SNI_output$stdN 	= NA*numeric(length(ages));
  SNI_output$CF 		= NA*numeric(length(ages));
  rawSNI 						= NA*numeric(length(ages));
  
  for(i in 1:length(ages)) { # Perform calculations at each sample age
    # Calculate window indexes
    if( i < round(0.5*(BandWidth/r))+1 ) {  # Samples near beginning (moving window truncated)
      SNI_output$winInd[i,] = c(1, ( i + round(0.5*(BandWidth/r)) )); 
    } else {
      if( i > length(ages)-round(0.5*(BandWidth/r)) ) { # Samples near end
        SNI_output$winInd[i,] = c( (i - round(0.5*(BandWidth/r))), length(ages) );  
      } else { 
        SNI_output$winInd[i,] = c( (i - round(0.5*(BandWidth/r))), (i+round(0.5*(BandWidth/r))) );
      }
    }
    
    # In each moving window...  
    
    # Get CHAR for samples in window (X is entire window population)
    X = CHAR[SNI_output$winInd[i,1]:SNI_output$winInd[i,2]]; 
    
    # inS, inN: boolean indexes of samples counted as S & N, respectively
    inS = (X >= thresh[SNI_output$winInd[i,1]:SNI_output$winInd[i,2]] );
    inN = !inS;
    
    # Fill in S & N populations, means, stds
    SNI_output$popS[[i]] = X[inS];
    SNI_output$popN[[i]] = X[inN];  
    
    if(is.na(mean(X[inN]))){
      SNI_output$meanN[i] = 0;  
    } else {
      SNI_output$meanN[i] = mean(X[inN]);
    }
    
    
    if(is.na(stats::sd(X[inN]))){
      SNI_output$stdN[i] = 1;
    } else {
      SNI_output$stdN[i] = stats::sd(X[inN]);
    }
    

    # Calculate correction factor
    v = length(SNI_output$popN[[i]]); # d.f. of N
    if(v == 0) {
      SNI_output$CF[i] = -1
    } else {
      SNI_output$CF[i] = (v - 2) / v; # SNI = Z * (v-2)/v  
    }
    
    
    # Calculate raw SNI (unsmoothed)
    if( length(SNI_output$popS[[i]])>0 ) {
      rawSNI[i] = (mean( (SNI_output$popS[[i]] - SNI_output$meanN[i]) ) 
                   / SNI_output$stdN[i]) * SNI_output$CF[i];
    } else {
      rawSNI[i] = 0; # SNI = 0 by definition when no samples exceed threshold
    }
    
  }
  
  # Smooth raw values to obtain final SNI
  SNI_output$SNI = stats::lowess(ages,rawSNI,f=(BandWidth/r)/length(ages),iter=0)$y;
  
  return(SNI_output);
  
}