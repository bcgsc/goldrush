import pandas as pd
from datetime import datetime
import sys
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from joblib import dump, load

def main():
    if len(sys.argv[1:]) != 2:
        print("Usage: {} <tsv> <file prefix>".format(sys.argv[0]))
        sys.exit()

    tsvfile = sys.argv[1]
    out_prefix = sys.argv[2]

    df = pd.read_csv(tsvfile, sep="\t")
    df = df.replace(to_replace="None", value=0)

    date_today = datetime.today()
    date_prefix = "{}-{}-{}".format(date_today.date(), date_today.month, date_today.year)

    df["length_over_yintercept"] = df["Length"]/df["yintercept"]

    df = df.replace(to_replace=float("inf"), value=0)

    X = df.loc[:, ["Correlation_coefficient", "yintercept", "slope", "num_mx",
                    "entropy", "mapped_bins", "length_over_yintercept"]].values

    # Feature Scaling
    sc = StandardScaler()
    X_train = sc.fit_transform(X)

    pca_analysis(X_train, df, date_prefix + "." + out_prefix)
    random_forest(X, df.loc[:, "Hairpin_status"].values)
    #tsne_analysis(X_train, df, date_prefix + "." + out_prefix)

def random_forest(X, y):
    "Try random forest classifier"
    from sklearn.model_selection import train_test_split
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    from sklearn.ensemble import RandomForestClassifier

    sc = load("scaler")
    #X_train = sc.fit_transform(X_train)
    #X_test = sc.transform(X_test)
    X = sc.transform(X)

    #classifier = RandomForestClassifier()
    #classifier.fit(X_train, y_train)
    classifier = load("random_forest_classifier")

    y_pred = classifier.predict(X)
    print(X)
    print(y_pred)

    from sklearn.metrics import confusion_matrix, accuracy_score, f1_score, precision_score
    print(confusion_matrix(y, y_pred))
    print(accuracy_score(y, y_pred))
    print(f1_score(y, y_pred, pos_label="Hairpin"))
    print(precision_score(y, y_pred, pos_label="Hairpin"))

    #dump(classifier, "random_forest_classifier")
    #dump(sc, "scaler")


def pca_analysis(X_train, df, prefix):
    # Applying PCA
    pca = PCA(n_components=2)
    X_train = pca.fit_transform(X_train)
    df1 = pd.DataFrame(X_train, columns=["PC1", "PC2"])
    colours = ["blue", "purple"]
    labels = ["Hairpin", "Non-hairpin"]
    for colour, label in zip(colours, labels):
        indices = df["Hairpin_status"] == label
        plt.scatter(x=df1.loc[indices, 'PC1'], y=df1.loc[indices, 'PC2'], c=colour, alpha=0.3)
    plt.legend(["Hairpin", "Non-hairpin"])
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA")

    plt.savefig("{}.pca_analysis.png".format(prefix), format="png")

def tsne_analysis(X_train, df, prefix):
    tsne = TSNE(init='pca')
    X_train = tsne.fit_transform(X_train)
    df1 = pd.DataFrame(X_train, columns=["PC1", "PC2"])

    colours = ["blue", "purple"]
    labels = ["Hairpin", "Non-hairpin"]
    for colour, label in zip(colours, labels):
        indices = df["Hairpin_status"] == label
        plt.scatter(x=df1.loc[indices, 'PC1'], y=df1.loc[indices, 'PC2'], c=colour, alpha=0.3)
    plt.legend(["Hairpin", "Non-hairpin"])
    plt.title("T-SNE")
    plt.xlabel("")
    plt.ylabel("")
    plt.savefig("{}.tsne_analysis.png".format(prefix), format="png")

if __name__ == "__main__":
    main()