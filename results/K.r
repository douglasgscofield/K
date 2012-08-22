reload <- function(doit=FALSE) if (doit) source("K.r")

form.filename <- function(U=1, s=1, h=0) 
{
    ans <- paste(sep="=", c("U", "s", "h"), c(U, s, h))
    ans <- paste(sep="_", "K", paste(collapse="_", ans))
    paste(sep="", ans, ".txt")
}

load <- function(U=1, s=1, h=0) 
{
    return(read.table(header=TRUE, sep="\t", file=form.filename(U=U, s=s, h=h)))
}

fig <- function(UU=c(0.02, 0.2, 1), s=1, h=0, do.postscript=FALSE)
{
    options(scipen=6)
    fig.lty <- seq(along=U)

    # load data
    df <- c()
    for (u in UU) {
        df <- rbind(df, load(U=u, s=s, h=h))
    }

    # set up plot
    #if (do.postscript) postscript(file="
    par(mfrow=c(2,2), pty="s", oma=c(2,1,1,1), mar=c(2,4,2,0), 
        mgp=c(2.5,0.3,0), las=1, tcl=0.3)

    # lethals
    plot(df$S_0, df$mean_totmuts, xlim=c(0,1), ylim=c(0.01, 1000), 
        xlab="", ylab="Mean lethals",
        type="n", log="y", frame.plot=TRUE)
    for (i in seq(along=UU)) {
        ddf <- subset(df, U == UU[i])
        lines(ddf$S_0, ddf$mean_totmuts, lty=fig.lty[i])
    }
    # inbreeding depression
    plot(df$S_0, df$IBD, xlim=c(0,1), ylim=c(0,1),
        xlab="", ylab="Inbreeding depression",
        type="n", log="", frame.plot=TRUE)
    for (i in seq(along=U)) {
        ddf <- subset(df, U == UU[i])
        lines(ddf$S_0, ddf$IBD, lty=fig.lty[i])
    }
    # secondary selfing rate
    plot(df$S_0, df$S_secondary, xlim=c(0,1), ylim=c(0,1),
        xlab="", ylab="Secondary selfing rate",
        type="n", log="", frame.plot=TRUE)
    for (i in seq(along=U)) {
        ddf <- subset(df, U == UU[i])
        lines(ddf$S_0, ddf$S_secondary, lty=fig.lty[i])
    }
    # Mean fitness
    plot(df$S_0, df$w_popmean, xlim=c(0,1), ylim=c(0,1),
        xlab="", ylab="Mean fitness",
        type="n", log="", frame.plot=TRUE)
    for (i in seq(along=U)) {
        ddf <- subset(df, U == UU[i])
        lines(ddf$S_0, ddf$w_popmean, lty=fig.lty[i])
    }

    mtext(paste(sep=":  ", "K output", 
        paste(collapse=" ... ", paste(sep="=", c("U","s","h"), 
        c(paste(collapse=", ", UU),s,h)))), outer=TRUE)
    mtext("Primary selfing rate", side=1, outer=TRUE)
}


