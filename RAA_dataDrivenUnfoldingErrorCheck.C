// Raghav Kunnawalkam Elayavalli
// Monday March 25th
// Rutgers

//
// Macro for data driven unfolding error correction since just plain unfolding doesnt fix the error bars correctly.
// here is the method:
//1) For each pT bin in your measured spectra, generate a gaussian distribution with the BinContent as mean and BinError as sigma. 
//2) generate a large number of spectra, say a thousand, by getting a random value from the gaussian for each pt bin. 
//3) now unfold these 1000 spectra, with the same response matrix (derived from recopt-genpt distributions)
//4) Now we have 1000 values of bin contents and bin errors for each pt bin. Fill the bin contents into a histogram and find the mean (which is our BinContent for our final unfolded histogram) and RMS (which the BinError) per pt bin. 
//
//

