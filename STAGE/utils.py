import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def setToArray(
        setInput,
        dtype='int64'
):
    """ This function transfer set to array.
        Args:
            setInput: set need to be trasnfered to array.
            dtype: data type.

        Return:
            arrayOutput: trasnfered array.
    """
    arrayOutput = np.zeros(len(setInput), dtype=dtype)
    index = 0
    for every in setInput:
        arrayOutput[index] = every
        index += 1
    return arrayOutput


def Euclidean_distance(
        feature_1,
        feature_2
):
    """ This function generates Euclidean distance between two vectors.
        Args:
            feature_1, feature_2: two vectors.

        Return:
            dist: distance between feature_1 and feature_2.
    """
    dist = np.linalg.norm(feature_1-feature_2)
    return dist


def dist_with_slice(
        new_spot,
        origin_coor_df
):
    """ This function generates Euclidean distance between two vectors.
        Args:
            new_spot: coordinate of spot.
            origin_coor_df: coordinate of all original spots.
        Return:
            min_dist_with_spot: minimun of distance between new spot and all original spots.
    """
    dist_with_spot=[]
    for it in range(origin_coor_df.shape[0]):
        dist_with_spot.append(Euclidean_distance(new_spot, origin_coor_df.iloc[it,:]))
    min_dist_with_spot=min(dist_with_spot)
    return min_dist_with_spot


def generation_coord_10x(
        adata,
        name='coord'
):
    """ This function generates spatial location for 10x Visium data.
        Args:
            adata: AnnData object storing original data.
            name: Item in adata.obsm used for get spatial location. Default is "coord".
        Return:
            coor_df: Spatial location of original data.
            fill_coor_df: Spatial location of generated data.
    """
    coor_df = pd.DataFrame(adata.obsm[name])

    coor_df_1 = pd.DataFrame(adata.obsm[name])
    coor_df_1.iloc[:, 1] = coor_df_1.iloc[:, 1] + 2 / 3

    coor_df_2 = pd.DataFrame(adata.obsm[name])
    coor_df_2.iloc[:, 0] = coor_df_2.iloc[:, 0] + 1
    coor_df_2.iloc[:, 1] = coor_df_2.iloc[:, 1] + 1 / 3

    coor_df_3 = pd.DataFrame(adata.obsm[name])
    coor_df_3.iloc[:, 0] = coor_df_3.iloc[:, 0] + 1
    coor_df_3.iloc[:, 1] = coor_df_3.iloc[:, 1] - 1 / 3

    fill_coor_df = pd.concat([coor_df, coor_df_1, coor_df_2, coor_df_3])
    fill_coor_df = fill_coor_df.drop_duplicates(subset=fill_coor_df.columns)

    coor_df.index=adata.obs.index
    coor_df.columns=["x","y"]
    fill_coor_df.columns = ["x", "y"]

    return coor_df, fill_coor_df


def generation_coord_ST(
        adata,
        name='coord'
):
    """ This function generates spatial location for Spatial Transcriptomics data.
        Args:
            adata: AnnData object storing original data.
            name: Item in adata.obsm used for get spatial location. Default is "coord".
        Return:
            coor_df: Spatial location of original data.
            fill_coor_df: Spatial location of generated data.
    """
    coor_df = pd.DataFrame(adata.obsm[name])

    coor_df_1 = pd.DataFrame(adata.obsm[name])
    coor_df_1.iloc[:, 1] = coor_df_1.iloc[:, 1] + 0.5

    coor_df_2 = pd.DataFrame(adata.obsm[name])
    coor_df_2.iloc[:, 0] = coor_df_2.iloc[:, 0] + 0.5

    coor_df_3 = pd.DataFrame(adata.obsm[name])
    coor_df_3.iloc[:, 0] = coor_df_3.iloc[:, 0] + 0.5
    coor_df_3.iloc[:, 1] = coor_df_3.iloc[:, 1] + 0.5

    fill_coor_df = pd.concat([coor_df, coor_df_1, coor_df_2, coor_df_3])
    fill_coor_df = fill_coor_df.drop_duplicates(subset=fill_coor_df.columns)

    coor_df.index=adata.obs.index
    coor_df.columns=["x","y"]
    fill_coor_df.columns = ["x", "y"]

    return coor_df, fill_coor_df


def recovery_coord(
        adata,
        name='coord',
        down_ratio=0.5
):
    """ This function generates spatial location for Spatial Transcriptomics data.
        Args:
            adata: AnnData object storing original data.
            name: Item in adata.obsm used for get spatial location. Default is "coord".
            down_ratio: Down-sampling ratio. Default is 0.5.
        Return:
            coor_df: Spatial location of dowm-sampled data.
            fill_coor_df: Spatial location of recovered data.
            sample_index: Index of downsampled data.
            sample_barcode: Barcode of downsampled data.
    """
    coor_df = pd.DataFrame(adata.obsm[name])
    coor_df.index = adata.obs.index
    coor_df.columns = ["x", "y"]
    sample_index=np.random.choice(range(coor_df.shape[0]), size=round(down_ratio*coor_df.shape[0]), replace=False)
    sample_index = setToArray(set(sample_index))
    sample_coor_df = coor_df.iloc[sample_index]
    sample_barcode = coor_df.index[sample_index]
    return sample_coor_df, coor_df, sample_index, sample_barcode


def Slide_seq_coord_3d(
        adata,
        name='coord',
        sec_name='section',
        select_section=[1, 3, 5, 6, 8],
        gap=0.05
):
    """ This function generates spatial location for Slide-seq data.
        Args:
            adata: AnnData object storing original data.
            name: Item in adata.obsm used for get spatial location. Default is "coord".
            sec_name: Item in adata.obs.columns used for choosing training sections.
            select_section: Index of training sections.
            gap: Distance between simulated and real sections. Half of distance between adjacent real sections.
        Return:
            coor_df: Spatial location of dowm-sampled data in real sections.
            fill_coor_df: Spatial location of recovered data in real sections.
            new_coor_df: Spatial location of generated data in simulated sections.
            all_coor_df: Spatial location of generated data in all sections.
            sample_xy_coord: 2d coordinates for simulated sections.
    """
    fill_coor_df = pd.DataFrame(adata.obsm[name])
    fill_coor_df.index = adata.obs.index
    fill_coor_df.columns = ['x', 'y', 'z']

    coor_df = fill_coor_df[adata.obs[sec_name].isin(select_section)]

    z_coord = fill_coor_df['z'].unique()
    new_z_coord = z_coord[range(z_coord.shape[0] - 1)] + gap

    border = round(np.max(np.amax(fill_coor_df[['x', 'y']], 0)) + 100)
    sample_xy_coord = pd.DataFrame(np.random.rand(border * 5, 2) * border, columns=["x", "y"])

    xmin, xmax = np.min(fill_coor_df['x']), np.max(fill_coor_df['x'])
    ymin, ymax = np.min(fill_coor_df['y']), np.max(fill_coor_df['y'])
    centre = [round(0.5 * (xmin + xmax)), round(0.5 * (ymin + ymax))]
    radius = min(round(0.5 * (ymax - ymin)), round(0.5 * (xmax - xmin)))

    dist_with_centre = sample_xy_coord[["x", "y"]].apply(lambda x: Euclidean_distance(x, centre), axis=1)
    sample_xy_coord = sample_xy_coord[dist_with_centre <= radius]

    new_coor_df = pd.DataFrame(columns=['x', 'y', 'z'])
    for i in new_z_coord:
        current_coor_df = pd.DataFrame(columns=['x', 'y'])
        current_coor_df[['x', 'y']] = sample_xy_coord[['x', 'y']]
        current_coor_df['z'] = i
        new_coor_df = pd.concat([new_coor_df, current_coor_df])

    all_coor_df = pd.concat([fill_coor_df, new_coor_df])

    return coor_df, fill_coor_df, new_coor_df, all_coor_df, sample_xy_coord


def show_train_hist(
        hist,
        loss_type,
        label,
        show = False,
        save = False,
        path = 'Train_hist.png'
):
    x = range(len(hist[loss_type]))

    y = hist[loss_type]

    plt.plot(x, y, label=label)

    plt.xlabel('Epoch')
    plt.ylabel('Loss')

    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout()

    if save:
        plt.savefig(path)

    if show:
        plt.show()
    else:
        plt.close()


