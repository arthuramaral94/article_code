# article_code
#ANN Training
#Using red and infrared (NIR) pixel values
#Normalized values
print(nntrad <- neuralnet(Z_Normalized ~ RED_N + NIR_N, data=Rtraditional,
hidden=3, rep=100, err.fct="sse", linear.output=TRUE))
plot(nntrad)
#Prediction for all other points
est_rnaREDNIR_TRAD_0901<- predict(nntrad,`dados`)
#Exporting the result of R to a .txt file
write.table(est_rnaREDNIR_TRAD_0901, file = "est_rnaREDNIR_TRAD_0901.txt",
row.names = TRUE, sep = ";", dec = ",")
