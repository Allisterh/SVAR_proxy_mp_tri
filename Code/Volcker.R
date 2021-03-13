Start = 1979.75
End1 = 1982.5
End2 = 1983.5

(T0A = Start - (End1 - Start) - 0.25)
(T0B = Start - (End2 - Start) - 0.25)
(T0C = T0A - (End1 - Start) - 0.25)
(TTA = End1 + (End1 - Start) + 0.25)
(TTB = End2 + (End2 - Start) + 0.25)


Volcker.fct <- function(x, series = 2, Par = NULL, Epsname = c("d", "s", "mp"), Start, End){
        if (!is.null(Par)){
       erg <- quasi.hd(x, series, Partial = Par, Epsname = Epsname, Start = Start, End = End, Freq = 4) 
        } else {
       erg <- c(quasi.hd(x, series, Partial = 1, Epsname = Epsname, Start = Start, End = End, Freq = 4),
                quasi.hd(x, series, Partial = 2, Epsname = Epsname, Start = Start, End = End, Freq = 4),
                quasi.hd(x, series, Partial = 3, Epsname = Epsname, Start = Start, End = End, Freq = 4)) %>% t %>% t                
        }
        erg
}



# Pi ----------------------------------------------------------------------
# CV 
pi_RefA_CV <- Volcker.fct(CV, series = 2, Start = Start, End = End1) #* 0.25
pi_RefB_CV <- Volcker.fct(CV, series = 2, Start = Start, End = End2) #* 0.25
pi_anteA_CV <- Volcker.fct(CV, series = 2, Start = T0A, End = Start- 0.25) #* 0.25
pi_anteB_CV <- Volcker.fct(CV, series = 2, Start = T0B, End = Start- 0.25) #* 0.25
pi_anteC_CV <- Volcker.fct(CV, series = 2, Start = T0C, End = T0A- 0.25) #* 0.25
pi_postA_CV <- Volcker.fct(CV, series = 2, Start = End1+ 0.25, End = TTA) #* 0.25
pi_postB_CV <- Volcker.fct(CV, series = 2, Start = End2+ 0.25, End = TTB) #* 0.25

# SW
pi_RefA_SW <- Volcker.fct(IV.SW, Par = 3, Start = Start, End = End1) #* 0.25
pi_RefB_SW <- Volcker.fct(IV.SW, Par = 3, Start = Start, End = End2) #* 0.25
pi_anteA_SW <- Volcker.fct(IV.SW, Par = 3, Start = T0A, End = Start- 0.25) #* 0.25
pi_anteB_SW <- Volcker.fct(IV.SW, Par = 3, Start = T0B, End = Start- 0.25) #* 0.25
pi_anteC_SW <- Volcker.fct(IV.SW, Par = 3, Start = T0C, End = T0A- 0.25) #* 0.25
pi_postA_SW <- Volcker.fct(IV.SW, Par = 3, Start = End1+ 0.25, End = TTA) #* 0.25
pi_postB_SW <- Volcker.fct(IV.SW, Par = 3, Start = End2+ 0.25, End = TTB) #* 0.25


# output ------------------------------------------------------------------
# CV
x_RefA_CV <- Volcker.fct(CV, series = 1, Start = Start, End = End1)  
x_RefB_CV <- Volcker.fct(CV, series = 1, Start = Start, End = End2)  
x_anteA_CV <- Volcker.fct(CV, series = 1, Start = T0A, End = Start- 0.25)  
x_anteB_CV <- Volcker.fct(CV, series = 1, Start = T0B, End = Start- 0.25)  
x_anteC_CV <- Volcker.fct(CV, series = 1, Start = T0C, End = T0A- 0.25)  
x_postA_CV <- Volcker.fct(CV, series = 1, Start = End1+ 0.25, End = TTA)  
x_postB_CV <- Volcker.fct(CV, series = 1, Start = End2+ 0.25, End = TTB)  

# SW
x_RefA_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = Start, End = End1)  
x_RefB_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = Start, End = End2)  
x_anteA_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = T0A, End = Start- 0.25)  
x_anteB_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = T0B, End = Start- 0.25)  
x_anteC_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = T0C, End = T0A- 0.25)  
x_postA_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = End1+ 0.25, End = TTA)  
x_postB_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = End2+ 0.25, End = TTB)  



