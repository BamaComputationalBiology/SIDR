import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.svm import NuSVC
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import model_selection as ms
from mpl_toolkits import mplot3d
import scipy.stats as sp
import os
if not os.path.exists('figures'):
    os.makedirs('figures')


target = []
not_target = []
unclassified = []


def import_seqdata(seqname, seqlen, gc, cov, hasBlast, isTarget, kdist_list):
    seqdata = (seqname, seqlen, kdist_list, isTarget, gc, cov)
    if hasBlast == 1:
        if isTarget == 1:
            target.append(seqdata)
        else:
            not_target.append(seqdata)
    else:
        unclassified.append(seqdata)
    return 0



def plot_mixtures(kmer_list):
    all_gc_t = []
    all_rcov_t = []
    all_kcov_t = []
    all_gc_n = []
    all_rcov_n = []
    all_kcov_n = []
    all_gc_u = []
    all_rcov_u = []
    all_kcov_u = []
    for i in range(len(kmer_list)):
        gc_t = []
        rcov_t = []
        kcov_t = []
        gc_n = []
        rcov_n = []
        kcov_n = []
        gc_u = []
        rcov_u = []
        kcov_u = []
        for seqdata in target:
            for j in range(len(seqdata[2][i])):
                if seqdata[2][i][j] == 0:
                    continue
                gc_t.append(seqdata[4])
                rcov_t.append(seqdata[5])
                kcov_t.append(j)
                all_gc_t.append(seqdata[4])
                all_rcov_t.append(seqdata[5])
                all_kcov_t.append(j)
        for seqdata in not_target:
            for j in range(len(seqdata[2][i])):
                if seqdata[2][i][j] == 0:
                    continue
                gc_n.append(seqdata[4])
                rcov_n.append(seqdata[5])
                kcov_n.append(j)
                all_gc_n.append(seqdata[4])
                all_rcov_n.append(seqdata[5])
                all_kcov_n.append(j)
        for seqdata in unclassified:
            for j in range(len(seqdata[2][i])):
                if seqdata[2][i][j] == 0:
                    continue
                gc_u.append(seqdata[4])
                rcov_u.append(seqdata[5])
                kcov_u.append(j)
                all_gc_u.append(seqdata[4])
                all_rcov_u.append(seqdata[5])
                all_kcov_u.append(j)
        plt.figure(i)
        ax = plt.axes(projection="3d")
        ax.scatter(gc_t, rcov_t, kcov_t, c='r', marker='o', alpha=0.1)
        ax.scatter(gc_u, rcov_u, kcov_u, c='g', marker='o', alpha=0.1)
        ax.scatter(gc_n, rcov_n, kcov_n, c='b', marker='o', alpha=0.1)
        plt.savefig(F"figures/mixture{kmer_list[i]}.png");
    plt.figure(len(kmer_list))
    ax = plt.axes(projection="3d")
    ax.scatter(all_gc_t, all_rcov_t, all_kcov_t, c='r', marker='o', alpha=0.1)
    ax.scatter(all_gc_u, all_rcov_u, all_kcov_u, c='g', marker='o', alpha=0.1)
    ax.scatter(all_gc_n, all_rcov_n, all_kcov_n, c='b', marker='o', alpha=0.1)
    plt.savefig("figures/mixtureALL.png");



def GM_model(kmer_list):
    totalHistData = []
    gm_models = {}
    aic_scores = {}
    bic_scores = {}
    for k in range(len(kmer_list)):
        kmerHistData = []
        for i in range(len(seq_list)):
            regionKmerHistData = readHistData(histFile_list[i][k])
            if blast_list[i] == 1 and tax_list[i] == 1:
                for dataPoint in regionKmerHistData:
                    kmerHistData.append(dataPoint)
            totalHistData.append((seq_list[i], kmer_list[k], regionKmerHistData))
        gm_models[kmer_list[k]] = [GaussianMixture(n, covariance_type='full').fit(kmerHistData) for n in range(1,6)]
        #aic_scores[kmer_list[k]] = [gm_models[kmer_list[k]][n].aic(kmerHistData) for n in range(0,5)]
        #bic_scores[kmer_list[k]] = [gm_models[kmer_list[k]][n].bic(kmerHistData) for n in range(0,5)]
    kmer_scores = {seq : {} for seq in seq_list}
    for seqTuple in totalHistData:
        kmer_scores[seqTuple[0]][seqTuple[1]] = [gm_models[seqTuple[1]][n].score(seqTuple[2]) for n in range(0,5)]
    return kmer_scores



def DT_classify(classifier):
    regionIDs = []
    X = []
    for item in unclassified:
        regionIDs.append(item[0])       # region name
        X.append(item[4:])
    Y = []
    for i in classifier.predict(X):
        Y.append(i)
    return list(zip(regionIDs, Y))



def DT_model():
    X = []
    Y = []
    features = ["GC", "Coverage"]
    for item in target + not_target:
        X.append(item[4:])
        Y.append(item[3])       # isTarget
    X_train, X_test, Y_train, Y_test = ms.train_test_split(X, Y, test_size=0.33, random_state=0)
    classifier = tree.DecisionTreeClassifier()
    classifier = classifier.fit(X_train, Y_train)
    print("Classifier built, score is %s out of 1.00" % classifier.score(X_test, Y_test))
    #with open("model.dot", 'w') as dotfile:
    #    tree.export_graphviz(classifier, out_file=dotfile, feature_names=features,
    #                         class_names=Y, filled=True, rounded=True, special_characters=True)
    return DT_classify(classifier)



def run_analysis(kmer_list):
    print("Seqdata imported, %d target sequences, %d contaminant sequences, and %d unclassified sequences\n" % (len(target), len(not_target), len(unclassified)))
    plot_mixtures(kmer_list)
## does this classifiers == even do anything????
    classifiers = [
##### here Case used a combination of paramaters from SVC and NuSVC. gamma goes to both, C only to SVC and nu only to NuSVC. so not sure what is correct. 
        NuSVC(gamma=2, nu=(len(target)/len(not_target))),
        #SVC(gamma=2, C=1, nu=(len(target)/len(not_target))),
        GaussianNB(),
        DecisionTreeClassifier(),
        MLPClassifier(alpha=1, max_iter=1000),
        RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        GaussianProcessClassifier(1.0 * RBF(1.0))]

## idk what test data is supposed to be. kinda seems like it's supposed to be all of the data since you would want to classify the data and then write it to the files? but not sure....
    if len(testdata) != 0:
        classifier = constructDTclassifier(corpus, kmer_list)
        result = classifyData(classifier, testdata)
    else:
        result = []
    with open("tokeep.txt", "w+") as k, open("toremove.txt", "w+") as r:
        for item in result:
            if item[1] == 1:
                k.write("%s x\n" % item[0])
            else:
                r.write("%s x\n" % item[0])
        for i in range(len(corpus)):
            if corpus[i][0] == 1:
                k.write("%s\n" % corpus_map[i])
            else:
                r.write("%s\n" % corpus_map[i])
    return 1
