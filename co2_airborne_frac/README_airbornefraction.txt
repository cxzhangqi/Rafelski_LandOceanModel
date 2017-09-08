README: detrendCO2inc2.m

This code finds the constant airborne fraction, and the airborne fraction anomaly.  

Some notes:
To get the airborne fraction anomaly of the model results, you'll need to save the results in a .mat file: save('filename.mat').  Then replace line 71 with load filename.mat.  I've included a sample .mat file (co2_update_VHMV_11jan11.mat) for you to test the code. 

Lines 97-104: The index numbers are hardcoded here (e.g. on line 101: ffi(1189:3104,2)), so you'll want to change the index numbers to work with your time range, or, to make your life easier later, determine the index numbers to use based on the "find" function (for examples, see lines 80-81 and 85-86) - then you can easily change the date range.  