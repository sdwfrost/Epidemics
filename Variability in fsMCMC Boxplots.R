#'
#'
#' Variability in fsMCMC Plots (N = 200, 1000)
#'


par(mfrow = c(1,1))

jpeg('beta_rep_vs_random.jpeg', width = 960, height = 960)
boxplot(as.numeric(Variability[1:100,3])*as.numeric(Variability[1:100, 2]),
        as.numeric(Variability[101:200,3])*as.numeric(Variability[101:200, 2]),
        as.numeric(Variability[201:300,3])*as.numeric(Variability[201:300, 2]),
        as.numeric(Variability[301:400,3])*as.numeric(Variability[301:400, 2]),
        ylab = expression(hat(psi)), names = c("Representative", "Random", "Representative", "Random"),
        ylim = c(0, 3.1), col = c("tomato2", "tomato2", "lightblue", "lightblue")
        )

legend("topright", legend = c(200, 1000), col = c("tomato2", "lightblue"), pch = 15,
       title = expression(N))

abline(a = 1, b = 0, col = "green", lty = 2)

dev.off()

jpeg('beta_rep_vs_random_1000.jpeg', width = 960, height = 960)
boxplot(
        as.numeric(Variability[201:300,3])*as.numeric(Variability[201:300, 2]),
        as.numeric(Variability[301:400,3])*as.numeric(Variability[301:400, 2]),
        ylab = expression(hat(psi)), names = c("Representative", "Random"),
        ylim = c(0.5, 1.5), col = c("lightblue", "lightblue")
)

legend("topright", legend = c(200, 1000), col = c("tomato2", "lightblue"), pch = 15,
       title = expression(N))

abline(a = 1, b = 0, col = "green", lty = 2)

dev.off()

jpeg("gamma_rep_vs_random.jpeg", width = 960, height = 960)
boxplot(as.numeric(Variability[1:100,5]),
        as.numeric(Variability[101:200,5]),
        as.numeric(Variability[201:300,5]),
        as.numeric(Variability[301:400,5]),
        ylab = expression(hat(gamma)), names = c("Representative", "Random", "Representative", "Random"),
        ylim = c(0, 0.4), col = c("tomato2", "tomato2", "lightblue", "lightblue")
        )


legend("topright", legend = c(200, 1000), col = c("tomato2", "lightblue"), pch = 15,
       title = expression(N))
abline(a = 0.15, b = 0, col = "red", lty = 2)
dev.off()


jpeg("gamma_rep_vs_random_1000.jpeg", width = 960, height = 960)
boxplot(
        as.numeric(Variability[201:300,5]),
        as.numeric(Variability[301:400,5]),
        ylab = expression(hat(gamma)), names = c("Representative", "Random"),
        ylim = c(0.1, 0.2), col = c("lightblue", "lightblue")
)


legend("topright", legend = c(200, 1000), col = c("tomato2", "lightblue"), pch = 15,
       title = expression(N))
abline(a = 0.15, b = 0, col = "red", lty = 2)
dev.off()


jpeg("ESS_200_vs_1000", width = 960, height = 960)
boxplot(as.numeric(Variability[1:100,9]),
        as.numeric(Variability[101:200,9]),
        as.numeric(Variability[201:300,9]),
        as.numeric(Variability[301:400,9]),
        ylab = "ESS", names = c("Representative", "Random", "Representative", "Random"),
        ylim = c(0, 450), col = c("tomato2", "tomato2", "lightblue", "lightblue")
        )

legend("topright", legend = c(200, 1000), col = c("tomato2", "lightblue"), pch = 15,
       title = expression(N))

dev.off()

jpeg("ESS_sec_200_vs_1000", width = 960, height = 960)
boxplot(as.numeric(Variability[1:100,9])/25,
        as.numeric(Variability[101:200,9])/25,
        as.numeric(Variability[201:300,9])/95,
        as.numeric(Variability[301:400,9])/95,
        ylim = c(0, 30),
        ylab = "ESS/sec", names = c("Representative", "Random", "Representative", "Random"),
        col = c("tomato2", "tomato2", "lightblue", "lightblue"))

legend("topright", legend = c(200, 1000), col = c("tomato2", "lightblue"), pch = 15,
       title = expression(N))

dev.off()
