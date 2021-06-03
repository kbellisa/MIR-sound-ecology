#####################################
#Global acoustic indices calculation#
#####################################

# 20160621

# Authors:   Amandine Gasc: 		amandine.gasc@gmail.com
#			Dante Francomano: 	dfrancomano@gmail.com



# spec<-matriceMoinsBruit_prop[,1]
# f=44100
# nmax = NULL
# amp =c(1/90,1/90)
# freq = NULL
# threshold = NULL
# plot = TRUE
# title = TRUE
# xlab = "Frequency (kHz)"
# ylab = "Amplitude"
# labels = TRUE
# legend = TRUE
# collab = "red" 



fpeaksFlat<-function (spec, f = NULL, nmax = NULL, amp = NULL, freq = 100, 
                      threshold = .30, plot = FALSE, title = FALSE, xlab = "Frequency (kHz)",
                      ylab = "Amplitude", labels = FALSE, legend = FALSE, collab = "red", 
                      ...) 
{
  if (is.matrix(spec)) {
    if (ncol(spec) != 2) 
      stop("If 'spec' is a numeric matrix it should be a two column matrix with the first colum describing the frequency x-axis and the second column describing the amplitude y-axis")
    N <- nrow(spec)
  }
  if (is.vector(spec)) {
    N <- length(spec)
    if (is.null(f)) {
      stop("If 'spec' is a numeric vector describing the amplitude only, the sampling frequency 'f' of the original signal should be provided (for instance (f = 44100)")
    }
    if (!is.null(f) && !is.na(f)) {
      spec <- cbind(seq(f/(N * 2), f/2, length = N)/1000, 
                    spec)
    }
    if (!is.null(f) && is.na(f)) {
      spec <- cbind(1:N, spec)
      plot <- FALSE
    }
  }
  flat <- round(N/20)
  spec.tmp <- c(spec[, 2], rep(NA, flat))
  for (i in 1:(N - (flat + 1))) {
    ref <- spec.tmp[i]
    for (j in 1:flat) {
      if (spec.tmp[i + j] == ref) {
        spec.tmp[i + j] <- spec.tmp[i + j] + 1e-05 * 
          spec.tmp[i + j]
      }
    }
  }
  spec <- cbind(spec[, 1], spec.tmp[1:N])
  sym <- discrets(spec[, 2], symb = 5, collapse = FALSE)
  
  # si ya des flats
  specF<-spec
  for (mod in 1:length(specF[,2]))
  {
    if(specF[mod,2]==0){specF[mod,2]<-specF[mod,2]+0.000001*mod}
  }
  sym2 <- discrets(specF[, 2], symb = 5, collapse = FALSE)
  #
  
  if (sym2[1] == "I") 
    sym2[1] <- "T"
  if (sym2[1] == "P") 
    sym2[1] <- "D"
  
  sym2 <- c(NA, sym2, NA)
  peaks <- which(sym2 == "P")
  valleys <- which(sym2 == "T")
  
  n <- length(peaks)
  if (n == 0) {
    res <- NA
    plot <- FALSE
  }    else {
    if (!is.null(amp) | !is.null(nmax)) {
      diffvp <- diffpv <- numeric(n)
      for (i in 1:n) {
        v <- specF[valleys[i], 2]
        p <- specF[peaks[i], 2]
        vv <- specF[valleys[i + 1], 2]
        diffvp[i] <- p - v
        diffpv[i] <- p - vv
      }
    }
    if (!is.null(nmax) && n != 0) {
      if (!is.null(amp) | !is.null(freq) | !is.null(threshold)) {
        cat("Caution! The argument 'nmax' overrides the arguments 'amp', 'freq', and 'threshold'")
      }
      if (n < nmax) {
        cat(paste("There are", n, "peaks only (< nmax ="), 
            nmax, ")")
      }
      if (nmax == 1) {
        tmp <- specF[peaks, , drop = FALSE]
        res <- tmp[which.max(tmp[, 2]), , drop = FALSE]
      }            else {
        alt <- cbind(peaks, diffvp, diffpv)
        leftorder <- alt[order(-alt[, 2]), , drop = FALSE]
        rightorder <- alt[order(-alt[, 3]), , drop = FALSE]
        left <- leftorder[, 1]
        right <- rightorder[, 1]
        l <- 0
        i <- 1
        while (l[i] < nmax) {
          comp <- left[1:i] %in% right[1:i]
          l <- c(l, length(comp[comp == TRUE]))
          i <- i + 1
        }
        peaks0 <- left[1:(i - 1)]
        if (l[i] > nmax) {
          error <- l[i] - nmax
          peaks0 <- peaks0[1:(length(peaks0) - error)]
        }
        peaks <- peaks0[comp]
        res <- matrix(na.omit(specF[peaks, ]), nc = 2)
        colnames(res) <- c("freq", "amp")
      }
    } else {
      if (!is.null(amp)) {
        if (length(amp) != 2) 
          stop("The length of 'amp' should equal to 2.")
        for (i in 1:n) {
          if (!is.na(diffvp[i]) && !is.na(diffpv[i]) && 
                diffvp[i] > 0 && diffpv[i] > 0 && diffvp[i] >= 
                amp[1] && diffpv[i] >= amp[2]) 
            peaks[i] <- peaks[i]
          else peaks[i] <- NA
        }
      }
      if (!is.null(freq)) {
        freq <- freq/1000
        diffpeak <- numeric(n - 1)
        for (i in 1:(n - 1)) {
          peak1 <- specF[peaks[i], 1]
          peak2 <- specF[peaks[i + 1], 1]
          diffpeak[i] <- peak2 - peak1
          if (!is.na(diffpeak[i]) && diffpeak[i] <= freq) 
            if (specF[peaks[i + 1], 2] > specF[peaks[i], 
                                               2]) 
            {peaks[i] <- NA}  else peaks[i + 1] <- NA
        }
      }
      if (!is.null(f) && is.na(f)) {
        res <- peaks
      }            else {
        res <- matrix(na.omit(specF[peaks, ]), nc = 2)
        colnames(res) <- c("freq", "amp")
      }
      if (!is.null(threshold)) {
        res <- res[res[, 2] > threshold, , drop = FALSE]
      }
    }
  }
  if (plot) {
    plot(spec, type = "l", xlab = xlab, ylab = ylab, xaxs = "i", 
         yaxt = "n", ...)
    if (title) {
      if (nrow(res) == 1) {
        text.title <- "peak detected"
      }            else {
        text.title <- "peaks detected"
      }
      title(main = paste(nrow(res), text.title))
    }
    points(res, col = collab)
    if (labels & nrow(res) != 0) 
      text(res, labels = round(res[, 1], 2), pos = 3, col = collab)
    if (!is.null(threshold)) {
      abline(h = threshold, col = collab, lty = 2)
      mtext(paste(threshold), side = 2, line = 0.5, at = threshold, 
            las = 1, col = collab)
    }
    if (legend) {
      if (!is.null(nmax)) {
        text.legend <- paste("nmax=", nmax, sep = "")
      }            else {
        if (is.null(amp)) {
          amp[1] <- amp[2] <- "-"
        }                else amp <- round(amp, 2)
        if (is.null(freq)) {
          freq <- "-"
        }
        if (is.null(threshold)) {
          threshold <- "-"
        }
        text.legend <- c(paste("amp=", amp[1], "/", amp[2], 
                               sep = ""), paste("freq=", freq, sep = ""), 
                         paste("threshold=", threshold, sep = ""))
      }
      legend("topright", pch = NA, legend = text.legend, 
             bty = "n", text.col = "darkgrey")
    }
    invisible(res)
  }    else return(res)
}



#############################
#############################     

CPSound <- function(wave,
                         leftChannel=TRUE,
                         rightChannel=TRUE,
                         min_freq = 0,
                         max_freq = 22050,
                         anthro_min = 0,
                         anthro_max = 1000,
                         bio_min=1000,
                         bio_max=22050,
                         wl=512,
                         j=5,
                         f=44100,
                         
                         AcouOccupancy=TRUE,
                         Bioac=TRUE,
                         Hf=TRUE,
                         Ht=TRUE,
                         H=FALSE,
                         ACI=TRUE,
                         BioLev_Joo=FALSE,
                         AEI_villa=TRUE,
                         M=TRUE,
                         AR=FALSE,
                         NDSI=TRUE,
                         ADI=TRUE,
                         NP=TRUE,
                         ASA=TRUE,
                         BN=TRUE,
                         SNR=TRUE,
                         AA=TRUE,
                         CAE=TRUE,
                         ADAE=TRUE,
                         Hm=FALSE,
                         Hv=FALSE,
                         MBA=FALSE,
                         spectradiv=FALSE,
                         spectralpersist=FALSE,
                         RMS=TRUE,
                         WIND=FALSE,
                         CRANE=FALSE,
                         CHORUS=FALSE)

{
  library(seewave)
  library(tuneR)
  library(soundecology)
  library(TTR)
  
  #arguments fpeak
  amp=c(1/90,1/90)
  freq=200
  plotpic=FALSE
  
  print("wave_before")
  str(wave)
  f<-wave@samp.rate
  print("wave_after")
  
  nyquist_freq <- f/2
  
  #Basic left channel calculations
  if(leftChannel==TRUE)
  {
    left <- mono(wave, which = c("left"))
    spec_left <- spectro(left, f = f, wl = wl, plot = FALSE, dB = "max0")$amp
    specA_left <- apply(spec_left, 1, meandB)
    rows_width = length(specA_left) / nyquist_freq
    min_row = round(min_freq * rows_width)
    max_row = round(max_freq * rows_width)
    specA_left_segment <- specA_left[min_row:max_row]
    freqs <- seq(from = min_freq, to = max_freq, length.out = length(specA_left_segment))
    specA_left_segment_normalized<-(specA_left_segment-min(specA_left_segment))/max(specA_left_segment-min(specA_left_segment))
    spec_L<-cbind(freqs,specA_left_segment_normalized)
  }
  
  #Basic right channel calculations
  if(rightChannel==TRUE)
  {
    right<- mono(wave, which = c("right"))
    spec_right <- spectro(right, f = f, wl = wl, plot = FALSE, dB = "max0")$amp
    specA_right <- apply(spec_right, 1, meandB)
    rows_width = length(specA_right) / nyquist_freq
    min_row = round(min_freq * rows_width)
    max_row = round(max_freq * rows_width)
    specA_right_segment <- specA_right[min_row:max_row]
    freqs <- seq(from = min_freq, to = max_freq, length.out = length(specA_right_segment))
    specA_right_segment_normalized<-(specA_right_segment-min(specA_right_segment))/max(specA_right_segment-min(specA_right_segment))
    spec_R<-cbind(freqs,specA_right_segment_normalized)
  }
  
  Table_left<-NULL
  Table_right<-NULL
  
  ###################################################################		Could not find ref ####### NEED TO CHANGE WEIGHTING !!??
  # # Accoustic Occupancy Index ref???
  if(AcouOccupancy==TRUE)
  {
    if(leftChannel==TRUE)
    {
      Spectrogram_Aleft<-spectro(left, f = f, wl = wl, plot = FALSE, dB = "A")$amp
      Spectrogram_Aleft_segment <- Spectrogram_Aleft[min_row:max_row, ]
      vec_left<-NULL
      for (i in 1:ncol(Spectrogram_Aleft_segment))
      {
        if(any(na.omit(Spectrogram_Aleft_segment[,i])>-30)){vec_left<-c(vec_left,1)} else{vec_left<-c(vec_left,0)}
      }
      AcouOccupancy_left<-length(vec_left[vec_left==1])/length(vec_left)
      Table_left<-cbind(Table_left,AcouOccupancy_left)
    }
    if(rightChannel==TRUE)
    {
      Spectrogram_Aright<-spectro(right, f = f, wl = wl, plot = FALSE, dB = "A")$amp
      Spectrogram_Aright_segment <- Spectrogram_Aright[min_row:max_row, ]
      vec_right<-NULL
      for (i in 1:ncol(Spectrogram_Aright_segment))
      {
        if(any(na.omit(Spectrogram_Aright_segment[,i])>-30)){vec_right<-c(vec_right,1)} else{vec_right<-c(vec_right,0)}
      }
      AcouOccupancy_right<-length(vec_right[vec_right==1])/length(vec_right)
      Table_right<-cbind(Table_right,AcouOccupancy_right)
    } 
  }
  
  ###################################################################
  # Relative avian abondance - (Boelman et al. 2007)
  #calculation of the dB en fonction to the frequency=meanspectrum db using
  #then calculation of the area under the spectrum dB
  if(Bioac==TRUE)
  {
    bioacouMeasure<-bioacoustic_index(wave, min_freq = min_freq, max_freq = max_freq, fft_w = wl)
    if(leftChannel==TRUE)
    {
      Bioac_left<-bioacouMeasure$left_area
      Table_left<-cbind(Table_left,Bioac_left)
    }
    if(rightChannel==TRUE)
    {
      Bioac_right<-bioacouMeasure$right_area
      Table_right<-cbind(Table_right,Bioac_right)
    }
  }
  
  
  ###################################################################
  #Temporal entropy - Ht (Sueur et al. 2008)
  if(Ht==TRUE)
  {
    th2<-function (env) 
    {
      options(warn = -1)
      if (is.na(env)) 
        z <- 0
      else 
      {
        options(warn = 0)
        if (any(env < 0, na.rm = TRUE)) 
          stop("data must be an envelope, i. e. a vector including positive values only.")
        N <- length(env)
        env <- env/sum(env)
        if (any(is.nan(env))) 
        {
          warning("Caution! There is no signal in this data set! The temporal entropy is null!", call. = FALSE)
          return(0)
        }
        if (sum(env)/(N * env[1]) == 1 | sum(env)/N == 1) 
        {
          warning("Caution! This is a square signal. The temporal entropy is null!", call. = FALSE)
          return(0)
        }
        env[env == 0] <- 1e-07
        env <- env/sum(env)
        z <- -sum(env * log(env))/log(N)
      }
      return(z)
    }
    if(leftChannel==TRUE)
    {
      env_left<-env(wave@left,f,plot=FALSE,envt="abs")
      env_left<-env_left/sum(as.numeric(env_left))
      Ht_left<-th2(env_left)
      Table_left<-cbind(Table_left,Ht_left)
    }
    if(rightChannel==TRUE)
    {
      env_right<-env(wave@right,f,plot=FALSE,envt="abs")
      env_right<-env_right/sum(as.numeric(env_right))
      Ht_right<-th2(env_right)
      Table_right<-cbind(Table_right,Ht_right)
    }
  }
  
  ###################################################################
  # Spectral entropy - Hf  (Sueur et al. 2008)
  if(Hf==TRUE)
  {
    if(leftChannel==TRUE)
    {
      Hf_left<-sh(spec_L)
      Table_left<-cbind(Table_left,Hf_left)
    }
    if(rightChannel==TRUE)
    {
      Hf_right<-sh(spec_R)
      Table_right<-cbind(Table_right,Hf_right)
    }
  }
  
  ###################################################################
  # Acoustic Entropy index - H (Sueur et al. 2008, joo et al. 2011)
  if(H==TRUE)
  {
    if(leftChannel==TRUE)
    {
      H_left=Hf_left*Ht_left
      Table_left<-cbind(Table_left,H_left)
    }
    if(rightChannel==TRUE)
    {
      H_right=Hf_right*Ht_right
      Table_right<-cbind(Table_right,H_right)
    }
  }
  
  ###################################################################		Could not find ref: Look for it in Sueur et al. 2014
  # Ratio of biophony to anthrophony - rho (Qi et al. 2008)
  # B=
  # A=
  # rho=B/A
  
  ###################################################################
  # ACI (Pieretti et al. 2011) ACI() {seewave} A REVOIR A PAROPOS DES FILTRES
  # Issue with the j with short sound# maybe by default =1 ???
  if(ACI==TRUE)
  {
    ACI_measure<-acoustic_complexity(wave,max_freq = max_freq, j = j, fft_w = wl)
    if(leftChannel==TRUE)
    {
      ACI_left<-ACI_measure$AciTotAll_left
      Table_left<-cbind(Table_left,ACI_left)
    }
    if(rightChannel==TRUE)
    {
      ACI_right<-ACI_measure$AciTotAll_right
      Table_right<-cbind(Table_right,ACI_right)
    }
  }
  
  ###################################################################     PLEASE CHECK
  # Bio and anthro level 1 (Joo et al., 2011)
  # vPSD_biophony #2 a 8 kHz
  # vPSD_anthrophony #1 a 2 kHz
  if(BioLev_Joo==TRUE)
  {
    if(leftChannel==TRUE)
    {
      Binspec_left<-soundscapespec(left,wl=wl,plot=FALSE)[,-1]
      Anthrolev_left<-Binspec_left[2]
      Biolev_left<-max(Binspec_left[3:8])
      Table_left<-cbind(Table_left,Anthrolev_left,Biolev_left)
    }
    if(rightChannel==TRUE)
    {
      Binspec_right<-soundscapespec(right,wl=wl,plot=FALSE)[,-1]
      Anthrolev_right<-Binspec_right[2]
      Biolev_right<-max(Binspec_right[3:8])
      Table_right<-cbind(Table_right,Anthrolev_right,Biolev_right)		
    }
  }
  
  ###################################################################		??? same as bio part of Joo et al. 2011?
  # Biophony level 2 (Krause et al. 2011)
  #bioPeak<-
  
  ###################################################################
  # Acoustic Evenness Index (AEI)== H'(Villanueva-Rivera et al. 2011) acoustic_evenness(TropicalSound)
  if(AEI_villa==TRUE)
  {
    AEI_villa_measure<-acoustic_evenness(wave, max_freq = max_freq, db_threshold = -36, freq_step = 1000)
    if(leftChannel==TRUE)
    {
      AEI_villa_left<-AEI_villa_measure$aei_left
      Table_left<-cbind(Table_left,AEI_villa_left)
    }
    if(rightChannel==TRUE)
    {
      AEI_villa_right<-AEI_villa_measure$aei_right
      Table_right<-cbind(Table_right,AEI_villa_right)
    }
  }
  
  ###################################################################     ??? Seems like it should be simpler than this
  # Median of amplitude envelope - M (Depraetere et al. 2012)
  if(M==TRUE)
  {
    bit<-wave@bit
    if(leftChannel==TRUE)
    {
      env_left<-env(wave@left,f,plot=FALSE,envt="abs")
      env_left<- env_left/sum(as.numeric(env_left)) 
      M_left<-median(env_left)*(2^(1-bit))
      Table_left<-cbind(Table_left,M_left)
    }
    if(rightChannel==TRUE)
    {
      env_right<-env(wave@right,f,plot=FALSE,envt="abs")
      env_right<- env_right/sum(as.numeric(env_right)) 
      M_right<-median(env_right)*(2^(1-bit))
      Table_right<-cbind(Table_right,M_right)
    }
  }
  
  ###################################################################		??? M is confusing and I don't get how we take the rank of a single value
  # Acoustic Richness - AR (Depraetere et al. 2012)
  if(AR==TRUE)
  {
    bit<-wave@bit
    th3<-function (env) 
    {
      options(warn = -1)
      if (is.na(env)) 
        z <- 0
      else 
      {
        options(warn = 0)
        if (any(env < 0, na.rm = TRUE)) 
          stop("data must be an envelope, i. e. a vector including positive values only.")
        N <- length(env)
        env <- env/sum(env)
        if (any(is.nan(env))) 
        {
          warning("Caution! There is no signal in this data set! The temporal entropy is null!", call. = FALSE)
          return(0)
        }
        if (sum(env)/(N * env[1]) == 1 | sum(env)/N == 1) 
        {
          warning("Caution! This is a square signal. The temporal entropy is null!", call. = FALSE)
          return(0)
        }
        env[env == 0] <- 1e-07
        env <- env/sum(env)
        z <- -sum(env * log(env))/log(N)
      }
      return(z)
    }
    if(leftChannel==TRUE)
    {
      env_left<-env(wave@left,f,plot=FALSE,envt="abs")
      env_left<- env_left/sum(as.numeric(env_left)) 
      Ht_left<-th3(env_left)
      M_left<-median(env_left)*(2^(1-bit))		
      AR_left<-(rank(M_left)*rank(Ht_left))/((duration(left)*left@samp.rate)^2)
      Table_left<-cbind(Table_left,AR_left)
    }
    if(rightChannel==TRUE)
    {
      env_right<-env(wave@right,f,plot=FALSE,envt="abs")
      env_right<- env_right/sum(as.numeric(env_right)) 
      Ht_right<-th3(env_right)
      M_right<-median(env_right)*(2^(1-bit))
      AR_right<-(rank(M_right)*rank(Ht_right))/((duration(right)*right@samp.rate)^2)
      Table_right<-cbind(Table_right,AR_right)
    }
  }
  # AR<-(rank(M)*rank(Ht))/(99^2)
  # Table<-cbind(Table,AR)
  # }
  
  ###################################################################     
  # Normalised difference soundscape index - NDSI (Kasten et al. 2012)
  if(NDSI==TRUE)
  {
    NDSI_measure<-ndsi(wave, fft_w = wl, anthro_min = anthro_min, anthro_max = anthro_max,bio_min = bio_min, bio_max = bio_max)
    if(leftChannel==TRUE)
    {
      NDSI_left<-NDSI_measure$ndsi_left
      #NDSI_seewave<-NDSI(soundscapespec(wave))
      Table_left<-cbind(Table_left,NDSI_left)
    }
    if(rightChannel==TRUE)
    {
      NDSI_right<-NDSI_measure$ndsi_right
      Table_right<-cbind(Table_right,NDSI_right)
    }
  }
  
  ###################################################################
  # Acoustic diversity index - ADI (=H') (Pekin et al. 2013)// H' (Villanueva-Riviera et al. 2011)
  if(ADI==TRUE)
  {
    ADI_measure<-acoustic_diversity(wave,max_freq = max_freq, db_threshold = -36, freq_step = 1000,shannon = TRUE)
    if(leftChannel==TRUE)
    {
      ADI_left<-ADI_measure$adi_left
      Table_left<-cbind(Table_left,ADI_left)
    }
    if(rightChannel==TRUE)
    {
      ADI_right<-ADI_measure$adi_right
      Table_right<-cbind(Table_right,ADI_right)
    }
  }
  
  ###################################################################
  # Sound pressure level parameters - L (Rychtarikova and Vermeir 2013)
  # This is not really an index more a parametrs and not linled with acoustic diversity...
  # Their parameters are all based on human perception of sound, making them anthropocentrically biased (A-weighting, critical bands, etc.).
  # But they do 2 cool things:
  # 1. Calculated max values of parameters for 1 - 99% of the time (using a 125 ms window for L) and plotted these for each recording
  # 2. Used maxes at 5, 50, and 95% for each parameter to cluster into "types" of soundscapes
  
  ###################################################################
  # Number of peaks - NP (Gasc et al. 2013b)
  if(NP==TRUE)
  {
    if(leftChannel==TRUE)
    {
      res1_left <- fpeaksFlat(spec_L,plot=FALSE,f)
      pictot_left<-nrow(res1_left)
      if(is.null(pictot_left)==FALSE)
      {
        if(pictot_left!=1)
        {
          res2_left<-fpeaksFlat(spec_L,amp=amp,freq=freq,plot=FALSE,f)
          npic_left<-nrow(res2_left)
        }
        else
        {
          res2_left<-fpeaksFlat(spec_L,plot=FALSE,f) #BG changed from spec_left
          npic_left<-nrow(res2_left)
        }
      }
      else
      {
        npic_left<-0
      }
      Table_left<-cbind(Table_left,npic_left)
    }
    if(rightChannel==TRUE)
    {
      res1_right <- fpeaksFlat(spec_R,plot=FALSE,f)
      pictot_right<-nrow(res1_right)
      if(is.null(pictot_right)==FALSE)
      {
        if(pictot_right!=1)
        {
          res2_right<-fpeaksFlat(spec_R,amp=amp,freq=freq,plot=FALSE,f)
          npic_right<-nrow(res2_right)
        }
        else
        {
          res2_right<-fpeaksFlat(spec_R,plot=FALSE,f)
          npic_right<-nrow(res2_right)
        }
      }
      else
      {
        npic_right<-0
      }
      Table_right<-cbind(Table_right,npic_right)
    }
  }
  
  ###################################################################
  # Average Signal Amplitude - ASA (Towsey et al., 2014)
  # normalised between 0 and 1
  # normalizing min/max = -50 dB and -3 dB
  if(ASA==TRUE)
  {
    if(leftChannel==TRUE)
    {
      dB_left<-(20*log10(mean(env(left,msmooth=c(512,0),plot=FALSE))/(2^(left@bit-1))))
      ASA_left<-((dB_left-(-36))/((-3)-(-36)))
      if(ASA_left<0)
      {
        ASA_left<-0
      }
      if(ASA_left>1)
      {
        ASA_left<-1
      }
      Table_left<-cbind(Table_left,ASA_left)
    }
    if(rightChannel==TRUE)
    {
      dB_right<-(20*log10(mean(env(right,msmooth=c(512,0),plot=FALSE))/(2^(right@bit-1))))
      ASA_right<-((dB_right-(-36))/((-3)-(-36)))
      if(ASA_right<0)
      {
        ASA_right<-0
      }
      if(ASA_right>1)
      {
        ASA_right<-1
      }
      Table_right<-cbind(Table_right,ASA_right)
    }	
  }
  
  ###################################################################
  # Background Noise - BN (Towsey et al., 2014; Towsey 2013)
  # normalised between 0 and 1
  # normalizing min/max = -50 dB and -3 dB
  if(BN==TRUE)
  {
    if(leftChannel==TRUE)
    {-30
      dB_left<-(20*log10(env(left,msmooth=c(512,0),plot=FALSE)/(2^(left@bit-1))))
      minenv_left<-min(dB_left)
      if(minenv_left<(-30))
      {
        minenv_left<-(-30)
        dB_left[dB_left<(-30)]<-(-30)
      }
      dBnoise_left<-dB_left[dB_left>=minenv_left&dB_left<=(minenv_left+10)]
      if(mean(dBnoise_left)!=(-30))
      {
        dBhist_left<-hist(dBnoise_left,breaks=99,plot=FALSE)
        dBvect_left<-as.vector(dBhist_left)
        smoothdB_left<-NULL
        smoothdB_left[1]<-dBvect_left$counts[1]
        smoothdB_left[100]<-dBvect_left$counts[100]
        for(k in 2:99)
        {
          smoothdB_left[k]<-((dBvect_left$counts[k-1]+dBvect_left$counts[k]+dBvect_left$counts[k+1])/3)
        }
        MI_left<-which.max(smoothdB_left)
        totalUnderMode_left<-sum(smoothdB_left[1:(MI_left-1)])
        SDthreshSum_left<-(0.68*totalUnderMode_left)
        SDcountSum_left<-NULL
        SDthresh_left<-NULL
        for(l in (MI_left-1):1)
        {
          SDcountSum_left<-sum(SDcountSum_left,smoothdB_left[l])
          if(SDcountSum_left>SDthreshSum_left)break
          SDthresh_left<-l
        }
        SD_left<-(MI_left-SDthresh_left)
        SDdB_left<-(0.1*SD_left)
        MIdB_left<-(minenv_left+(0.1*MI_left))
        N_left<-0
        BNdB_left<-(MIdB_left+(N_left*SDdB_left))
      }else{
        BNdB_left<-(-30)
      }
      BN_left<-((BNdB_left-(-30))/((-3)-(-30)))
      if(BN_left<0)
      {
        BN_left<-0
      }
      if(BN_left>1)
      {
        BN_left<-1
      }
      Table_left<-cbind(Table_left,BN_left)
    }
    if(rightChannel==TRUE)
    {
      dB_right<-(20*log10(env(right,msmooth=c(512,0),plot=FALSE)/(2^(right@bit-1))))
      minenv_right<-min(dB_right)
      if(minenv_right<(-30))
      {
        minenv_right<-(-30)
        dB_right[dB_right<(-30)]<-(-30)
      }
      dBnoise_right<-dB_right[dB_right>=minenv_right&dB_right<=(minenv_right+10)]
      if(mean(dBnoise_right)!=(-30))
      {
        dBhist_right<-hist(dBnoise_right,breaks=99,plot=FALSE)
        dBvect_right<-as.vector(dBhist_right)
        smoothdB_right<-NULL
        smoothdB_right[1]<-dBvect_right$counts[1]
        smoothdB_right[100]<-dBvect_right$counts[100]
        for(k in 2:99)
        {
          smoothdB_right[k]<-((dBvect_right$counts[k-1]+dBvect_right$counts[k]+dBvect_right$counts[k+1])/3)
        }
        MI_right<-which.max(smoothdB_right)
        totalUnderMode_right<-sum(smoothdB_right[1:(MI_right-1)])
        SDthreshSum_right<-(0.68*totalUnderMode_right)
        SDcountSum_right<-NULL
        SDthresh_right<-NULL
        for(l in (MI_right-1):1)
        {
          SDcountSum_right<-sum(SDcountSum_right,smoothdB_right[l])
          if(SDcountSum_right>SDthreshSum_right)break
          SDthresh_right<-l
        }
        SD_right<-(MI_right-SDthresh_right)
        SDdB_right<-(0.1*SD_right)
        MIdB_right<-(minenv_right+(0.1*MI_right))
        N_right<-0
        BNdB_right<-(MIdB_right+(N_right*SDdB_right))
      }else{
        BNdB_right<-(-30)
      }
      BN_right<-((BNdB_right-(-30))/((-3)-(-30)))
      if(BN_right<0)
      {
        BN_right<-0
      }
      if(BN_right>1)
      {
        BN_right<-1
      }
      Table_right<-cbind(Table_right,BN_right)
    }
  }

  ###################################################################
  # Signal-to-Noise Ratio - (Towsey et al., 2014)
  # normalised between 0 and 1
  # normalizing min/max = 3 dB and 50 dB
  if(SNR==TRUE)
  {
    if(leftChannel==TRUE)
    {
      dB_left<-(20*log10(env(left,msmooth=c(512,0),plot=FALSE)/(2^(left@bit-1))))
      maxenv_left<-max(dB_left)
      minenv_left<-min(dB_left)
      if(minenv_left<(-36))
      {
        minenv_left<-(-36)
        dB_left[dB_left<(-36)]<-(-36)
      }
      dBnoise_left<-dB_left[dB_left>=minenv_left&dB_left<=(minenv_left+10)]
      if(mean(dBnoise_left)!=(-36))
      {
        dBhist_left<-hist(dBnoise_left,breaks=99,plot=FALSE)
        dBvect_left<-as.vector(dBhist_left)
        smoothdB_left<-NULL
        smoothdB_left[1]<-dBvect_left$counts[1]
        smoothdB_left[100]<-dBvect_left$counts[100]
        for(k in 2:99)
        {
          smoothdB_left[k]<-((dBvect_left$counts[k-1]+dBvect_left$counts[k]+dBvect_left$counts[k+1])/3)
        }
        MI_left<-which.max(smoothdB_left)
        totalUnderMode_left<-sum(smoothdB_left[1:(MI_left-1)])
        SDthreshSum_left<-(0.68*totalUnderMode_left)
        SDcountSum_left<-NULL
        SDthresh_left<-NULL
        for(l in (MI_left-1):1)
        {
          SDcountSum_left<-sum(SDcountSum_left,smoothdB_left[l])
          if(SDcountSum_left>SDthreshSum_left)break
          SDthresh_left<-l
        }
        SD_left<-(MI_left-SDthresh_left)
        SDdB_left<-(0.1*SD_left)
        MIdB_left<-(minenv_left+(0.1*MI_left))
        N_left<-0
        BNdB_left<-(MIdB_left+(N_left*SDdB_left))
      }else{
        BNdB_left<-(-36)
      }
      SNRdB_left<-(maxenv_left-BNdB_left)
      SNR_left<-((SNRdB_left-3)/(36-3))
      if(SNR_left<0)
      {
        SNR_left<-0
      }
      if(SNR_left>1)
      {
        SNR_left<-1
      }
      Table_left<-cbind(Table_left,SNR_left)
    }
    if(rightChannel==TRUE)
    {
      dB_right<-(20*log10(env(right,msmooth=c(512,0),plot=FALSE)/(2^(right@bit-1))))
      maxenv_right<-max(dB_right)
      minenv_right<-min(dB_right)
      if(minenv_right<(-36))
      {
        minenv_right<-(-36)
        dB_right[dB_right<(-36)]<-(-36)
      }
      dBnoise_right<-dB_right[dB_right>=minenv_right&dB_right<=(minenv_right+10)]
      if(mean(dBnoise_right)!=(-36))
      {
        dBhist_right<-hist(dBnoise_right,breaks=99,plot=FALSE)
        dBvect_right<-as.vector(dBhist_right)
        smoothdB_right<-NULL
        smoothdB_right[1]<-dBvect_right$counts[1]
        smoothdB_right[100]<-dBvect_right$counts[100]
        for(k in 2:99)
        {
          smoothdB_right[k]<-((dBvect_right$counts[k-1]+dBvect_right$counts[k]+dBvect_right$counts[k+1])/3)
        }
        MI_right<-which.max(smoothdB_right)
        totalUnderMode_right<-sum(smoothdB_right[1:(MI_right-1)])
        SDthreshSum_right<-(0.68*totalUnderMode_right)
        SDcountSum_right<-NULL
        SDthresh_right<-NULL
        for(l in (MI_right-1):1)
        {
          SDcountSum_right<-sum(SDcountSum_right,smoothdB_right[l])
          if(SDcountSum_right>SDthreshSum_right)break
          SDthresh_right<-l
        }
        SD_right<-(MI_right-SDthresh_right)
        SDdB_right<-(0.1*SD_right)
        MIdB_right<-(minenv_right+(0.1*MI_right))
        N_right<-0
        BNdB_right<-(MIdB_right+(N_right*SDdB_right))
      }else{
        BNdB_right<-(-36)
      }
      SNRdB_right<-(maxenv_right-BNdB_right)
      SNR_right<-((SNRdB_right-3)/(36-3))
      if(SNR_right<0)
      {
        SNR_right<-0
      }
      if(SNR_right>1)
      {
        SNR_right<-1
      }
      Table_right<-cbind(Table_right,SNR_right)
    }
  }
  
  ###################################################################
  # Acoustic Activity - (Towsey et al., 2014)
  # normalised between 0 and 1
  if(AA==TRUE)
  {
    if(leftChannel==TRUE)
    {
      dB_left<-(20*log10(env(left,msmooth=c(512,0),plot=FALSE)/(2^(left@bit-1))))
      minenv_left<-min(dB_left)
      if(minenv_left<(-36))
      {
        minenv_left<-(-36)
        dB_left[dB_left<(-36)]<-(-36)
      }
      dBnoise_left<-dB_left[dB_left>=minenv_left&dB_left<=(minenv_left+10)]
      if(mean(dBnoise_left)!=(-36))
      {
        dBhist_left<-hist(dBnoise_left,breaks=99,plot=FALSE)
        dBvect_left<-as.vector(dBhist_left)
        smoothdB_left<-NULL
        smoothdB_left[1]<-dBvect_left$counts[1]
        smoothdB_left[100]<-dBvect_left$counts[100]
        for(k in 2:99)
        {
          smoothdB_left[k]<-((dBvect_left$counts[k-1]+dBvect_left$counts[k]+dBvect_left$counts[k+1])/3)
        }
        MI_left<-which.max(smoothdB_left)
        totalUnderMode_left<-sum(smoothdB_left[1:(MI_left-1)])
        SDthreshSum_left<-(0.68*totalUnderMode_left)
        SDcountSum_left<-NULL
        SDthresh_left<-NULL
        for(l in (MI_left-1):1)
        {
          SDcountSum_left<-sum(SDcountSum_left,smoothdB_left[l])
          if(SDcountSum_left>SDthreshSum_left)break
          SDthresh_left<-l
        }
        SD_left<-(MI_left-SDthresh_left)
        SDdB_left<-(0.1*SD_left)
        MIdB_left<-(minenv_left+(0.1*MI_left))
        N_left<-0
        BNdB_left<-(MIdB_left+(N_left*SDdB_left))
      }else{
        BNdB_left<-(-36)
      }
      SNRperFrame_left<-(dB_left-BNdB_left)
      AA_left<-((length(which(SNRperFrame_left>3)))/length(dB_left))
      
      Table_left<-cbind(Table_left,AA_left)
    }
    if(rightChannel==TRUE)
    {
      dB_right<-(20*log10(env(right,msmooth=c(512,0),plot=FALSE)/(2^(right@bit-1))))
      minenv_right<-min(dB_right)
      if(minenv_right<(-36))
      {
        minenv_right<-(-36)
        dB_right[dB_right<(-36)]<-(-36)
      }
      dBnoise_right<-dB_right[dB_right>=minenv_right&dB_right<=(minenv_right+10)]
      if(mean(dBnoise_right)!=(-36))
      {
        dBhist_right<-hist(dBnoise_right,breaks=99,plot=FALSE)
        dBvect_right<-as.vector(dBhist_right)
        smoothdB_right<-NULL
        smoothdB_right[1]<-dBvect_right$counts[1]
        smoothdB_right[100]<-dBvect_right$counts[100]
        for(k in 2:99)
        {
          smoothdB_right[k]<-((dBvect_right$counts[k-1]+dBvect_right$counts[k]+dBvect_right$counts[k+1])/3)
        }
        MI_right<-which.max(smoothdB_right)
        totalUnderMode_right<-sum(smoothdB_right[1:(MI_right-1)])
        SDthreshSum_right<-(0.68*totalUnderMode_right)
        SDcountSum_right<-NULL
        SDthresh_right<-NULL
        for(l in (MI_right-1):1)
        {
          SDcountSum_right<-sum(SDcountSum_right,smoothdB_right[l])
          if(SDcountSum_right>SDthreshSum_right)break
          SDthresh_right<-l
        }
        SD_right<-(MI_right-SDthresh_right)
        SDdB_right<-(0.1*SD_right)
        MIdB_right<-(minenv_right+(0.1*MI_right))
        N_right<-0
        BNdB_right<-(MIdB_right+(N_right*SDdB_right))
      }else{
        BNdB_right<-(-36)
      }
      SNRperFrame_right<-(dB_right-BNdB_right)
      AA_right<-((length(which(SNRperFrame_right>3)))/length(dB_right))
      
      Table_right<-cbind(Table_right,AA_right)
    }
  }
  
  ###################################################################
  # Count of Acoustic Events - (Towsey et al., 2014)
  # normalised between 0 and 1
  # normalizing min/max = 0 and 140
  if(CAE==TRUE)
  {
    if(leftChannel==TRUE)
    {
      dB_left<-(20*log10(env(left,msmooth=c(512,0),plot=FALSE)/(2^(left@bit-1))))
      minenv_left<-min(dB_left)
      if(minenv_left<(-30))
      {
        minenv_left<-(-30)
        dB_left[dB_left<(-30)]<-(-30)
      }
      dBnoise_left<-dB_left[dB_left>=minenv_left&dB_left<=(minenv_left+10)]
      if(mean(dBnoise_left)!=(-30))
      {
        dBhist_left<-hist(dBnoise_left,breaks=99,plot=FALSE)
        dBvect_left<-as.vector(dBhist_left)
        smoothdB_left<-NULL
        smoothdB_left[1]<-dBvect_left$counts[1]
        smoothdB_left[100]<-dBvect_left$counts[100]
        for(k in 2:99)
        {
          smoothdB_left[k]<-((dBvect_left$counts[k-1]+dBvect_left$counts[k]+dBvect_left$counts[k+1])/3)
        }
        MI_left<-which.max(smoothdB_left)
        totalUnderMode_left<-sum(smoothdB_left[1:(MI_left-1)])
        SDthreshSum_left<-(0.68*totalUnderMode_left)
        SDcountSum_left<-NULL
        SDthresh_left<-NULL
        for(l in (MI_left-1):1)
        {
          SDcountSum_left<-sum(SDcountSum_left,smoothdB_left[l])
          if(SDcountSum_left>SDthreshSum_left)break
          SDthresh_left<-l
        }
        SD_left<-(MI_left-SDthresh_left)
        SDdB_left<-(0.1*SD_left)
        MIdB_left<-(minenv_left+(0.1*MI_left))
        N_left<-0
        BNdB_left<-(MIdB_left+(N_left*SDdB_left))
      }else{
        BNdB_left<-(-30)
      }
      SNRperFrame_left<-(dB_left-BNdB_left)
      CAEraw_left<-numeric(length=1)
      NumberFrames_left<-(length(SNRperFrame_left))
      for(m in (3:(NumberFrames_left+1)))
      {
        if((((m-2)==1)&(SNRperFrame_left[1]>3))|(SNRperFrame_left[m-1]>3&SNRperFrame_left[m-2]<=3))
        {
          CAEraw_left<-sum(CAEraw_left,1)
        }
      }
      CAE_left<-((CAEraw_left-0)/(140-0))
      if(CAE_left<0)
      {
        CAE_left<-0
      }
      if(CAE_left>1)
      {
        CAE_left<-1
      }
      Table_left<-cbind(Table_left,CAE_left)
    }
    if(rightChannel==TRUE)
    {
      dB_right<-(20*log10(env(right,msmooth=c(512,0),plot=FALSE)/(2^(right@bit-1))))
      minenv_right<-min(dB_right)
      if(minenv_right<(-30))
      {
        minenv_right<-(-30)
        dB_right[dB_right<(-30)]<-(-30)
      }
      dBnoise_right<-dB_right[dB_right>=minenv_right&dB_right<=(minenv_right+10)]
      if(mean(dBnoise_right)!=(-30))
      {
        dBhist_right<-hist(dBnoise_right,breaks=99,plot=FALSE)
        dBvect_right<-as.vector(dBhist_right)
        smoothdB_right<-NULL
        smoothdB_right[1]<-dBvect_right$counts[1]
        smoothdB_right[100]<-dBvect_right$counts[100]
        for(k in 2:99)
        {
          smoothdB_right[k]<-((dBvect_right$counts[k-1]+dBvect_right$counts[k]+dBvect_right$counts[k+1])/3)
        }
        MI_right<-which.max(smoothdB_right)
        totalUnderMode_right<-sum(smoothdB_right[1:(MI_right-1)])
        SDthreshSum_right<-(0.68*totalUnderMode_right)
        SDcountSum_right<-NULL
        SDthresh_right<-NULL
        for(l in (MI_right-1):1)
        {
          SDcountSum_right<-sum(SDcountSum_right,smoothdB_right[l])
          if(SDcountSum_right>SDthreshSum_right)break
          SDthresh_right<-l
        }
        SD_right<-(MI_right-SDthresh_right)
        SDdB_right<-(0.1*SD_right)
        MIdB_right<-(minenv_right+(0.1*MI_right))
        N_right<-0
        BNdB_right<-(MIdB_right+(N_right*SDdB_right))
      }else{
        BNdB_right<-(-30)
      }
      SNRperFrame_right<-(dB_right-BNdB_right)
      CAEraw_right<-numeric(length=1)
      NumberFrames_right<-(length(SNRperFrame_right))
      for(m in (3:(NumberFrames_right+1)))
      {
        if((((m-2)==1)&(SNRperFrame_right[1]>3))|(SNRperFrame_right[m-1]>3&SNRperFrame_right[m-2]<=3))
        {
          CAEraw_right<-sum(CAEraw_right,1)
        }
      }
      CAE_right<-((CAEraw_right-0)/(140-0))
      if(CAE_right<0)
      {
        CAE_right<-0
      }
      if(CAE_right>1)
      {
        CAE_right<-1
      }
      Table_right<-cbind(Table_right,CAE_right)
    }
  }
  
  ###################################################################
  # Average Duration of Acoustic Events - (Towsey et al., 2014)
  # normalised between 0 and 1
  # normalizing min/max = 0 and 500 milliseconds
  if(ADAE==TRUE)
  {
    if(leftChannel==TRUE)
    {
      dB_left<-(20*log10(env(left,msmooth=c(512,0),plot=FALSE)/(2^(left@bit-1))))
      minenv_left<-min(dB_left)
      if(minenv_left<(-30))
      {
        minenv_left<-(-30)
        dB_left[dB_left<(-30)]<-(-30)
      }
      dBnoise_left<-dB_left[dB_left>=minenv_left&dB_left<=(minenv_left+10)]
      if(mean(dBnoise_left)!=(-30))
      {
        dBhist_left<-hist(dBnoise_left,breaks=99,plot=FALSE)
        dBvect_left<-as.vector(dBhist_left)
        smoothdB_left<-NULL
        smoothdB_left[1]<-dBvect_left$counts[1]
        smoothdB_left[100]<-dBvect_left$counts[100]
        for(k in 2:99)
        {
          smoothdB_left[k]<-((dBvect_left$counts[k-1]+dBvect_left$counts[k]+dBvect_left$counts[k+1])/3)
        }
        MI_left<-which.max(smoothdB_left)
        totalUnderMode_left<-sum(smoothdB_left[1:(MI_left-1)])
        SDthreshSum_left<-(0.68*totalUnderMode_left)
        SDcountSum_left<-NULL
        SDthresh_left<-NULL
        for(l in (MI_left-1):1)
        {
          SDcountSum_left<-sum(SDcountSum_left,smoothdB_left[l])
          if(SDcountSum_left>SDthreshSum_left)break
          SDthresh_left<-l
        }
        SD_left<-(MI_left-SDthresh_left)
        SDdB_left<-(0.1*SD_left)
        MIdB_left<-(minenv_left+(0.1*MI_left))
        N_left<-0
        BNdB_left<-(MIdB_left+(N_left*SDdB_left))
      }else{
        BNdB_left<-(-30)
      }
      SNRperFrame_left<-(dB_left-BNdB_left)
      CAEraw_left<-numeric(length=1)
      NumberFrames_left<-(length(SNRperFrame_left))
      for(m in (3:(NumberFrames_left+1)))
      {
        if((((m-2)==1)&(SNRperFrame_left[1]>3))|(SNRperFrame_left[m-1]>3&SNRperFrame_left[m-2]<=3))
        {
          CAEraw_left<-sum(CAEraw_left,1)
        }
      }
      AEframes_left<-sum(SNRperFrame_left>3[TRUE])
      if(CAEraw_left!=0)
      {
        ADAEraw_left<-((1000*512/left@samp.rate)*(AEframes_left/CAEraw_left))
        ADAE_left<-((ADAEraw_left-0)/(500-0))
        if(ADAE_left<0)
        {
          ADAE_left<-0
        }
        if(ADAE_left>1)
        {
          ADAE_left<-1
        }
      }else{
        ADAE_left<-0
      }
      Table_left<-cbind(Table_left,ADAE_left)
    }
    if(rightChannel==TRUE)
    {
      dB_right<-(20*log10(env(right,msmooth=c(512,0),plot=FALSE)/(2^(right@bit-1))))
      minenv_right<-min(dB_right)
      if(minenv_right<(-30))
      {
        minenv_right<-(-30)
        dB_right[dB_right<(-30)]<-(-30)
      }
      dBnoise_right<-dB_right[dB_right>=minenv_right&dB_right<=(minenv_right+10)]
      if(mean(dBnoise_right)!=(-30))
      {
        dBhist_right<-hist(dBnoise_right,breaks=99,plot=FALSE)
        dBvect_right<-as.vector(dBhist_right)
        smoothdB_right<-NULL
        smoothdB_right[1]<-dBvect_right$counts[1]
        smoothdB_right[100]<-dBvect_right$counts[100]
        for(k in 2:99)
        {
          smoothdB_right[k]<-((dBvect_right$counts[k-1]+dBvect_right$counts[k]+dBvect_right$counts[k+1])/3)
        }
        MI_right<-which.max(smoothdB_right)
        totalUnderMode_right<-sum(smoothdB_right[1:(MI_right-1)])
        SDthreshSum_right<-(0.68*totalUnderMode_right)
        SDcountSum_right<-NULL
        SDthresh_right<-NULL
        for(l in (MI_right-1):1)
        {
          SDcountSum_right<-sum(SDcountSum_right,smoothdB_right[l])
          if(SDcountSum_right>SDthreshSum_right)break
          SDthresh_right<-l
        }
        SD_right<-(MI_right-SDthresh_right)
        SDdB_right<-(0.1*SD_right)
        MIdB_right<-(minenv_right+(0.1*MI_right))
        N_right<-0
        BNdB_right<-(MIdB_right+(N_right*SDdB_right))
      }else{
        BNdB_right<-(-30)
      }
      SNRperFrame_right<-(dB_right-BNdB_right)
      CAEraw_right<-numeric(length=1)
      NumberFrames_right<-(length(SNRperFrame_right))
      for(m in (3:(NumberFrames_right+1)))
      {
        if((((m-2)==1)&(SNRperFrame_right[1]>3))|(SNRperFrame_right[m-1]>3&SNRperFrame_right[m-2]<=3))
        {
          CAEraw_right<-sum(CAEraw_right,1)
        }
      }
      AEframes_right<-sum(SNRperFrame_right>3[TRUE])
      if(CAEraw_right!=0)
      {
        ADAEraw_right<-((1000*512/right@samp.rate)*(AEframes_right/CAEraw_right))
        ADAE_right<-((ADAEraw_right-0)/(500-0))
        if(ADAE_right<0)
        {
          ADAE_right<-0
        }
        if(ADAE_right>1)
        {
          ADAE_right<-1
        }
      }else{
        ADAE_right<-0
      }
      Table_right<-cbind(Table_right,ADAE_right)
    }	
  }
  
  ###################################################################
  # Mid-band activity - (Towsey et al., 2014) the fraction of spectrogram cells in themid-band (482 Hz–3500 Hz) where the spectral amplitude exceeds 0.015.
  if(MBA==TRUE)
  {
    if(leftChannel==TRUE)
    {
      spectrogram_left<-spectro(left,dB=NULL)
      freqBinThresh_left<-NULL
      for (i in(1:length(spectrogram_left$freq)))
      {
        specRowi_left<-spectrogram_left$amp[i,]
        roundLengthSpecRowi<-round_any((length(specRowi_left)/8),8)
        specHist_left<-hist(specRowi_left,breaks=seq(min(min(specRowi_left)),max(max(specRowi_left)),length=roundLengthSpecRowi), plot = FALSE)
        histVect_left<-as.vector(specHist_left)
        smoothHist_left<-NULL
        smoothHist_left[1]<-histVect_left$counts[1]
        smoothHist_left[2]<-((histVect_left$counts[1]+histVect_left$counts[2]+histVect_left$counts[3])/3)
        smoothHist_left[roundLengthSpecRowi-2]<-((histVect_left$counts[roundLengthSpecRowi-2]+histVect_left$counts[roundLengthSpecRowi-3]+histVect_left$counts[roundLengthSpecRowi-4])/3)
        smoothHist_left[roundLengthSpecRowi-1]<-histVect_left$counts[roundLengthSpecRowi-1]
        for(k in 3:(roundLengthSpecRowi-3))
        {
          smoothHist_left[k]<-((histVect_left$counts[k-2]+histVect_left$counts[k-1]+histVect_left$counts[k]+histVect_left$counts[k+1]+histVect_left$counts[k+1])/5)
        }
        if(which.max(smoothHist_left)/(roundLengthSpecRowi-1)>0.95)
        {
          MI_left<-floor(0.95*(roundLengthSpecRowi-1))
        }else{
          MI_left<-which.max(smoothHist_left)
        }
        MIdB_left<-specHist_left$mids[MI_left]
        totalUnderMode_left<-sum(smoothHist_left[1:(MI_left-1)])
        SDthreshSum_left<-(0.68*totalUnderMode_left)
        SDcountSum_left<-NULL
        SDthresh_left<-0
        if(MI_left>1)
        {
          for(l in (MI_left-1):1)
          {
            SDcountSum_left<-sum(SDcountSum_left,smoothHist_left[l])
            if(SDcountSum_left>SDthreshSum_left)break
            SDthresh_left<-l
            SD_left<-(MIdB_left-histVect_left$mids[SDthresh_left])
          }
        }else{
          SDthresh_left<-min(smoothHist_left)
          SD_left<-(MIdB_left-histVect_left$mids[1])
        }
        N<-0.05
        freqBinThresh_left[i]<-(MIdB_left+(N*SD_left))
      }
      smoothFreqBinThresh_left<-NULL
      smoothFreqBinThresh_left[1]<-freqBinThresh_left[1]
      smoothFreqBinThresh_left[2]<-((freqBinThresh_left[1]+freqBinThresh_left[2]+freqBinThresh_left[3])/3)
      smoothFreqBinThresh_left[length(freqBinThresh_left)-1]<-((freqBinThresh_left[length(freqBinThresh_left)-2]+freqBinThresh_left[length(freqBinThresh_left)-1]+freqBinThresh_left[length(freqBinThresh_left)])/3)
      smoothFreqBinThresh_left[length(freqBinThresh_left)]<-freqBinThresh_left[length(freqBinThresh_left)]
      for(m in 3:(length(freqBinThresh_left)-2))
      {
        smoothFreqBinThresh_left[m]<-((freqBinThresh_left[m-2]+freqBinThresh_left[m-1]+freqBinThresh_left[m]+freqBinThresh_left[m+1]+freqBinThresh_left[m+1])/5)
      }
      
      spectrogramRed1_left<-spectrogram_left
      for(n in(1:length(spectrogram_left$freq)))
      {
        spectrogramRed1_left$amp[n,]<-(spectrogram_left$amp[n,]-smoothFreqBinThresh_left[n])
      }
    }
    

    midBandAmp_left<-spectro(left)$amp[spectro(left)$freq<=3.500&spectro(left)$freq>=0.65,]
    MBA_left<-(length(midBandAmp_left[midBandAmp_left>0.015])/length(midBandAmp_left))
    
    SPECTRO<-spectro(wave)$amp[spectro(wave)$freq<=3.500&spectro(wave)$freq>=0.65,]
    SPECTRO2<-length(SPECTRO>0.015)/length(SPECTRO)
  }
  
  ###################################################################
  # Entropy of spectral maxima - Hm (Towsey et al., 2014)
  # if(Hm==TRUE)
  # {
  # source("J:/programmes/R/analyse son/R-seewave/code fonction generalisee/all acoustic index/Hm.r")
  # Hm_index<-Hm(wave,f)
  # Table<-cbind(Table,Hm)
  # }
  
  # ###################################################################
  # # Entropy of spectral variance - Hv (Towsey et al., 2014)
  # #a la palce d un spectre de moyenne on a un spectre de variance
  # if(Hv==TRUE)
  # {
  # source("F:/projets/En_cours/Arizona_2014_Purdue/acoustic measure/Hv.r")
  # Hv<-Hv(wave)
  # Table<-cbind(Table,Hv)
  # }
  
  # ###################################################################
  # # Spectral diversity - (Towsey et al., 2014)
  # if(spectradiv==TRUE)
  # {
  # source("F:/projets/En_cours/Arizona_2014_Purdue/acoustic measure/sptradiv.r")
  # spectradiv<-spectradiv(wave,f)
  # Table<-cbind(Table,spectradiv)
  # }
  
  # ###################################################################
  # # Spectral persistence - (Towsey et al., 2014)
  # if(spectralpersist==TRUE)
  # {
  # source("F:/projets/En_cours/Arizona_2014_Purdue/acoustic measure/spectralpersist.r")
  # spectralpersist<-spectralpersist(wave,f)
  # Table<-cbind(Table,spectralpersist)
  # }
  
  ###################################################################
  # # RMS (and RMS L-R difference)
  
  if(RMS==TRUE)
  {
    if(leftChannel==TRUE)
    {
      RMS_left<-((rms(oscillo(left, plot = FALSE))))
      Table_left<-cbind(Table_left,RMS_left)
    }
    
    if(rightChannel==TRUE)
    {
      RMS_right<-((rms(oscillo(right, plot = FALSE))))
      Table_right<-cbind(Table_right,RMS_right)
    }
    if(leftChannel==TRUE & rightChannel==TRUE)
    {
      rmsLminusR<-(RMS_left-RMS_right)
      Table_left<-cbind(Table_left,rmsLminusR)
      Table_right<-cbind(Table_right,rmsLminusR)
    }
  }
  

############# WIND ANALYSIS ###########

if(WIND==TRUE)
{
  if(leftChannel==TRUE)
  {
    mean.left <- meanspec(left,f=f,wl = wl,plot=FALSE, norm = FALSE,PSD=TRUE)
    cut.left <- cutspec(mean.left,flim=c(0,.5))
    cut.left <- colMeans(cut.left)
    WIND_left <- cut.left[2]
    Table_left<-cbind(Table_left,WIND_left)
  }
  
  if(rightChannel==TRUE)

  {
    mean.right <- meanspec(right,f=f,wl = wl,plot=FALSE, norm = FALSE,PSD = TRUE)
    cut.right <- cutspec(mean.right,flim=c(0,.5))
    cut.right <- colMeans(cut.right)
    WIND_right <- cut.right[2]
    Table_right <-cbind(Table_right,WIND_right)
  }
  }

############# CRANE ANALYSIS ###########

if(CRANE==TRUE)
{
  if(leftChannel==TRUE)
  {
    cut.left <- cutspec(mean.left,flim=c(.65,1.3))
    cut.left <- colMeans(cut.left)
    CRANE_left <- cut.left[2]
    Table_left<-cbind(Table_left,CRANE_left)
  }
  
  if(rightChannel==TRUE)
    
  {
    cut.right <- cutspec(mean.right,flim=c(.65,1.3))
    cut.right <- colMeans(cut.right)
    CRANE_right <- cut.right[2]
    Table_right <-cbind(Table_right,CRANE_right)
  }
}

############# CHORUS ANALYSIS ###########

if(CHORUS==TRUE)
{
  if(leftChannel==TRUE)
  {
    cut.left <- cutspec(mean.left,flim=c(2.7,3.4))
    cut.left <- colMeans(cut.left)
    CHORUS_LEFT <- cut.left[2]
    Table_left<-cbind(Table_left,CHORUS_LEFT)
  }
  
  if(rightChannel==TRUE)
    
  {
    cut.right <- cutspec(mean.right,flim=c(2.7,3.4))
    cut.right <- colMeans(cut.right)
    CHORUS_right <- cut.right[2]
    Table_right <-cbind(Table_right,CHORUS_right)
  }
}


########## PREPARE DATA FRAME ############ 

Table_left<-as.data.frame(Table_left)
Table_right<-as.data.frame(Table_right)

print(Table_left)
print(Table_right)

if(leftChannel==FALSE)
{
  Table_left[]<-NA
}
if(rightChannel==FALSE)   
{
  Table_right[]<-NA
}
Result<-list(Table_right)
names(Result)<-c("Mono_right")
return(Result)
}