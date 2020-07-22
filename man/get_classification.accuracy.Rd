\name{get.classification.accuracy}
\alias{get.classification.accuracy}
\title{get.classification.accuracy}
\description{
This function performs classification accuracy analysis using the training and 
test sets. Users can choose support vector machine, logistic regression,
random forest, and naive bayes as classifiers. The performance evaluation
is determined based on the total classification rate, balanced accuracy rate,
and AUC.  
}
\usage{
get.classification.accuracy(kfold, featuretable, classlabels, 
kernelname = "radial", errortype = "AUC", conflevel = 95,
 classifier = "svm", seednum = 555,
 testfeaturetable = NA, testclasslabels = NA)
}
\arguments{
  \item{kfold}{
  	Number of folds for cross-validation. e.g. 5 or 10
  }
  \item{featuretable}{
R object for feature table with only differentially expressed features. 
This is the training set. The first two columns should be m/z and time.
}
  \item{classlabels}{
Class labels vector. e.g. c("case","control","case")
}
  \item{kernelname}{
Kernel for SVM: e.g. "radial" or "linear"
}
  \item{errortype}{
total: total classification accuracy rate;
(number of correct classifications/total N)
BAR: balanced accuracy rate; accounts of number of 
correct classification per class; 
BAR=(1/m)*((C1/N1)+(C2/N2)+...+(Cm/Nm))
where m is the number of classes, Cm is the number of
correct classifications in class m, and Nm is the total
number of subjects in class m.
AUC: area under the curve
}
  \item{conflevel}{
Confidence level for k-fold classification accuracy
e.g: 95
}
  \item{classifier}{
Classification algorithm to be used for ROC analysis. 
svm: Support Vector Machine
logitreg: Logistic Regression
rf: random forest
naivebayes: naive bayes
eg: "svm", "logitreg", "rf", "naivebayes"
}
  \item{seednum}{
Starting point used in the generation of a sequence of
random numbers. e.g. 555
}
  \item{testfeaturetable}{
R object for test feature table with only differentially expressed features. 
This is the test set. The first two columns should be m/z and time. The
order of features should be same as the training set.
}
  \item{testclasslabels}{
Class labels vector for samples in the test set.
}
}
\details{
Function to evaluate classification.
This function performs classification accuracy analysis using the training and 
test sets. Users can choose support vector machine, logistic regression,
random forest, and naive bayes as classifiers. The performance evaluation
is determined based on the total classification rate, balanced accuracy rate,
and AUC. 

}
\value{
Classification accuracy in training and test sets
}
\author{
Karan Uppal; kuppal2@emory.edu
}
