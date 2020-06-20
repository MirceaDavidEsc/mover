##### Discretization or Subsetting Functions ######

#' Discretize a trajectory by temporal steps
#' Note: this function is to be used if there are unequal steps between measurements. If the time step
#' is uniform you should use subsetData as it is much faster
#' @param track A trajectory data frame with the first column representing time.
#' @param threshT The time step to use when discretizing.
#'
#' @return A subset of trajectory points with appropriate temporal discretization.
#' @export
#'
#' @examples
discretizeTemporally = function(track, threshT) {
  track$keep = FALSE
  lastTime = track[1,1]
  track$keep[1] = T
  for (thisFrame in 1:dim(track)[1]) {
    if (track[thisFrame,1] >= lastTime + threshT) {
      track$keep[thisFrame] = T
      lastTime = track[thisFrame,1]

    }
  }
  discrTrack = track[track$keep,]
  discrTrack = discrTrack[,-ncol(discrTrack)]
  return(discrTrack)
}

###FUNCTION###
#INPUT:
# track: a data frame that contains at least three columns (with the correct order required): Frame/Time, X position, Y position
# threshL: The threshold length of discretization.
#OUTPUT: discretizedTrack is a subset of the input track, with only the row entries selected such that the distance between consecutive positions in entry rows equals or exceeds the threshold limit.
discretizeSpatially <- function(track,threshL) {
  discretizedTrack <- track[1,]
  initialPos <- track[1,2:3]
  for (currFrame in 1:length(track[,1]))
  {
    #print(paste("Frame",currFrame))
    nextPos <- track[currFrame,2:3]
    distance <- sqrt(sum((nextPos - initialPos)^2))
    #print(distance)
    if (distance > threshL)
    {
      #print("ADDED!")
      initialPos <- track[currFrame,2:3]
      discretizedTrack <- rbind(discretizedTrack,track[currFrame,])
    }
  }
  return(discretizedTrack)
}


#' Quick subset of data frames
#' Quickly subset a data frame to every nth entry, where n is defined by nskip.
#' @param dataFrame
#' @param nskip
#'
#' @return
#' @export
#'
#' @examples
subsetData = function(dataFrame, nskip) {
  indsToKeep = seq(from = 1, to = nrow(dataFrame), by = nskip)
  return(dataFrame[indsToKeep,])
}


##### Basic trajectory measurements #####

####FUNCTION:
#Receives a time series of x and y components of a vector. fromPrevious specifies if the function returns absolute angle
#From (1,0) vector, or if you want the angle relative to the previous vector.
getAngles <- function(x,y,fromPrevious = FALSE)
{
  angles <- c(rep(0,length(x)))
  #v <- normalize(x,y)

  #print(length(angles))
  if (!fromPrevious)
  {
    angles <-atan2(y,x)
  }
  else
  {
    baseAngle <- atan2(y[1],x[1])
    for(i in 1:length(x))
    {
      newAngle <- atan2(y[i],x[i])
      diff <- newAngle - baseAngle
      if (is.na(diff)) { # If the difference is NA, then just assign NA.
        angles[i] = newAngle
      } else { # need to calculate the new angle

        if (diff > pi) #Was some leftward-turn that went over the negative pi - pi line
        {
          diff <- (-2*pi+diff)
        } else if(diff < -pi) # Rightward turn crossing the -pi pi line
        {
          diff <- 2*pi+diff
        }
        angles[i] <- diff
      }

      baseAngle <- newAngle
    }
    angles = c(angles[2:length(angles)],NA)
  }
  return(angles)
}

####FUNCTION: From a given input of time points and positions, it calculates the speeds
calculateSpeed = function(trajectorydata)
{
  colnames(trajectoryData) = c("Frame", "X", "Y")
  trajectorydata$dX = c(NA, diff())
  speed = distance/displacement[,1]
  return(speed)
}



#' Measure point-by-point displacement in trajectory
#' Measures
#'
#' @param movementData the original data frame of positions, discreized appropriately.
#' @param colNames the names you wish to assign the columns for the displacements data frame.
#'
#' @return Returns a data frame with colNames column names with change in time and position.
#' @export
#'
#' @examples
calculateDisplacement = function(movementData,colNames=c('deltaT','deltaX','deltaY','deltaZ')) {
  numPoints = nrow(movementData)
  dimensionality = ncol(movementData)
  startingPoints = movementData[1:numPoints-1,]
  endingPoints = movementData[2:numPoints,]
  displacement = endingPoints - startingPoints
  colnames(displacement) = colNames[1:dimensionality]
  return(displacement)
}



#Calculates the displacements, speeds, absolute vector angles, and relative angles of movement data.
#INPUT:
#* movementData: a N-by-M dimensional vector containing, in each row, the time of recording and the coordinate values of the point at that time.
deriveTrajectoryData = function(movementData) {
  dim = length(movementData[1,])
  displacement = calculateDisplacement(movementData)
  absAngle = getAngles(displacement, fromPrevious = FALSE)
  absAngle = c(absAngle,NA)
  relAngle = getAngles(displacement, fromPrevious = TRUE)
  relAngle = c(relAngle,NA)
  displacement = rbind(displacement,c(rep(NA,length(displacement[1,]))))
  speed = calculateSpeed(displacement)

  movementData = cbind(movementData,displacement)
  movementData$Speed = speed
  movementData$AbsoluteAngle = absAngle
  movementData$RelativeAngle = relAngle
  return(movementData)
}


#' Calculate local tortuosity at any point in a trajectory
#'
#' @param trajectory A data frame containing the step-by-step displacements of a trajectory (dt, dX, dY)
#' @param bufferSize The frame range over which to measure tortuosity.
#'
#' @return
#' @export
#'
#' @examples
measureTortuosity = function(trajectory,bufferSize)
{
  if (bufferSize < 1 | bufferSize > nrow(trajectory)/2)  {
    print("Buffer size is invalid.")
    return
  }
  minIndex = bufferSize+1
  maxIndex = nrow(trajectory)-bufferSize

  tortuosity = c(rep(NA,nrow(trajectory)))

  steps = calculateDisplacement(trajectory)
  colnames(steps) = c("deltaX", "deltaY")
  stepDistance = c(0,sqrt((steps$deltaX^2 + steps$deltaY^2)))
  cumStepDistance = cumsum(stepDistance)

  for (i in minIndex:maxIndex)  {
    #Calculate distance between start and end frame.
    startPos = trajectory[i-bufferSize,]
    endPos = trajectory[i+bufferSize,]
    distance = sqrt(sum((startPos-endPos)^2))

    #Calculate the sum over the entire path length of the trajectory.
    totalDistance = cumStepDistance[i+bufferSize]-cumStepDistance[i-bufferSize]
    #Calculate the tortuosity and save it into
    tortuosity[i] = distance/totalDistance
  }
  return(tortuosity)
}

###FUNCTION: Takes in a time series of centroid positions for an animal and plots a
# square plot of the trajectory.
plotTrajectory = function(centerData)
{
  maxBound = max(max(abs(centerData$CentroidX)),max(abs(centerData$CentroidY)))
  trajectoriesPx = ggplot(centerData,aes(x=CentroidX,y=CentroidY)) + geom_point() + xlab("Centroid X Coordinate (px)") + ylab("Centroid Y Coordinate (px)") + coord_fixed(ratio=1,xlim = c(-maxBound,maxBound), ylim=c(-maxBound,maxBound)) + theme_bw()
  return(trajectoriesPx)
}


averageByBreaks <- function(dataFrame,breaks) {
  # This function takes a dataframe and a list of breaks, and creates averages of the entries between each break.
  meanVals = c(rep(0,length(breaks)))
  lastBreak = -Inf
  for (i in 2:length(breaks)) {
    meanVals[i-1] = mean(dataFrame[dataFrame[,1] > lastBreak & dataFrame[,1] <  breaks[i],2])
    lastBreak = breaks[i]
  }
  meanVals[length(breaks)] = NA
  return(meanVals)
}


##### Statistical trajectory measures #####
##FUNCTION
# INPUT: A series of turning angles that you wish to sequentially compare.
# OUTPUT: A contingency table of LL, LR, RL, and RR turn sequences, as well as total L and R turns.
anglePairContingencyTable <- function(turningAngles)
{
  turns <- c(rep(0,6)) ##[LL,LR,RL,RR,Lturns,Rturns]
  for (i in 2:length(turningAngles))
  {
    if (turningAngles[i] > 0) #Left turn
    {
      turns[5] <- turns[5]+1
      if (turningAngles[i-1] > 0) #Following a Left Turn
      {
        turns[1] <- turns[1]+1
      } else  { #Following a right turn
        turns[3] <- turns[3]+1
      }
    } else { #Right turn
      turns[6] <- turns[6]+1
      if (turningAngles[i-1] > 0) #following a left turn
      {
        turns[2] <- turns[2]+1
      } else #Following a right turn
      {
        turns[4] <- turns[4]+1
      }
    }
  }
  Total <- length(turningAngles)-1
  expectedLL <- (turns[5]*turns[5])/(Total^2)
  expectedRR <- (turns[6]*turns[6])/(Total^2)
  expectedLR <- (turns[5]*turns[6])/(Total^2)
  expectedRL <- (turns[6]*turns[5])/(Total^2)
  return(c(turns,expectedLL,expectedLR,expectedRL,expectedRR))
}


shortTimeACF = function(data,windowSize,offset,maxLag=windowSize/2)
{
  windowSeq = seq(1,length(data)-windowSize-1,offset)
  acfData = matrix(0,length(windowSeq),maxLag+1)

  for (i in 1:length(windowSeq))
  {
    start = windowSeq[i]
    subData = data[start:(start+windowSize)]
    subtMean = subData - mean(subData)
    autocorr = acf(subtMean,lag.max=maxLag,plot=FALSE)
    acfData[i,] = autocorr$acf
  }
  return(acfData)
}


##### Mean Square Displacement Functions #####

meanSquareDisplacement <- function(individualTrack,step,numRuns) {
  ###FUNCTION###
  #Calculates mean square displacement from a track data file.
  #INPUT:
  # individualTrack: a N-by-x data frame that has a time, X, and Y column in them.
  # step: the time-step to take between recorded MSD's.
  # numRuns: the number of runs over which to subset a track.

  #OUTPUT:
  # storeMSD: A matrix of mean square displacements for all time points chosen, from 1 to maxDelta. Also stores the frame number, standard deviation, and sample size per time point.


  #CALCULATING MEAN-SQUARE DISPLACEMENT BETWEEN RUNS
  #Find the longest list to set as the run length
  print("Calculating minimum and maximum lengths")

  storeMSDData(individualTrack,step,numRuns)

  uniqueRuns = unique(individualTrack$run)
  for (i in uniqueRuns)
  {
    print(i)
    run = as.factor(i)
    selectedTrack <- individualTrack[individualTrack$run==run,]
    frame <-individualTrack[individualTrack$run==run,1] #Extract all the relevant information about a particular run.
    X <- individualTrack[individualTrack$run==run,2]
    Y <- individualTrack[individualTrack$run==run,3]
    for (time in 1:length(timePoints))
    {
      #print(paste("time",time))
      finalPos <- c(selectedTrack$X_um[selectedTrack$Time_S==(timePoints[time])],selectedTrack$Y_um[selectedTrack$Time_S==(timePoints[time])])
      initialPos <- c(selectedTrack$X_um[selectedTrack$Time_S==0],selectedTrack$X_um[selectedTrack$Time_S==0])
      if (!length(finalPos) == 0)
      {
        displace <- sum(finalPos-initialPos)
        storeSD[time,run] <- sum(displace^2)
      }
    }
  }
  for (i in 1:length(timePoints)) #Get actual MSD data from the square displacement.
  {
    dataToUse <- na.omit(storeSD[i,])
    storeMSD$time[i]<-(timePoints[i])
    storeMSD$MSD[i] <- mean(dataToUse)
    storeMSD$StdDev[i] <- sd(dataToUse)
    storeMSD$numEntries[i] <- length(dataToUse)
  }
  storeMSD$StdDev[is.na(storeMSD$StdDev)] <- 0
  return(storeMSD)
}

StoreMSDData = function(individualTrack,step,numRuns)
{
  maxLength <- max(individualTrack[,1])
  timePoints = seq(0,maxLength,step) #The time steps ate which to record MSD
  numEntries <- length(timePoints) #How many data points are we taking?
  print(paste('numEntries',numEntries,''))
  storage <- c(rep(0,numEntries))
  storeMSD <- data.frame(time=storage,MSD=storage,StdDev=storage,numEntries=storage) #Create data frame to store entries
  storeSD <- matrix(NA,numEntries,numRuns)
}



###FUNCTION: Calculates the mean square displacement within an inidividual trajectory.It does so by dividing the trajectory into numRuns runs and calculating the MSD at the time points step distance apart.
#INPUT:
# individualTrack:
# step:
# numRuns:
#OUTPUT:
MSDWithinRun = function(individualTrack,step,numRuns)
{
  print("Calculating MSD by sub-sampling run.")
  runLength <- floor(max(individualTrack[,1])/numRuns) #The length of each track.

  startFrames <- seq(0,max(individualTrack[,1]),runLength) #The sequence of starting times when sub-dividing track into multiple independent runs.
  startFrames = startFrames[1:numRuns]

  timePoints = seq(0,runLength,step) #The series of dT times at which to measure square displacement within a run.
  numEntries <- length(timePoints) #The total sample size of runs created by the following track separation.
  storeSD <- matrix(0,numEntries,numRuns) #Storage for square displacement in each individual run.
  storage <- c(rep(0,numEntries))
  storeMSD <- data.frame(time=storage,MSD=storage,StdDev=storage,numEntries=storage) #Create data frame to store entries

  measureTimes = sapply(startFrames,function(data) newData=data+timePoints) %>% as.matrix()+1
  cxs = (individualTrack$CentroidX[measureTimes]) %>% matrix(nrow=numEntries,ncol=numRuns) %>% t()
  cxs2 = (cxs - individualTrack$CentroidX[measureTimes[1,]])^2
  cys = (individualTrack$CentroidY[measureTimes]) %>% matrix(nrow=numEntries,ncol=numRuns) %>% t()
  cys2 = (cys - individualTrack$CentroidY[measureTimes[1,]])^2

  SDStore = cxs2+cys2
  MSD = apply(X = SDStore,MARGIN = 2,mean)

  #print(paste("recordFrame",recordFrame,"Time_s",individualTrack$Time_S))
  storeMSD$time<-timePoints
  storeMSD$MSD <- apply(X = SDStore,MARGIN = 2,mean)
  storeMSD$StdDev <- apply(X = SDStore,MARGIN = 2,sd)
  storeMSD$numEntries <- numRuns
  return(storeMSD)
}
