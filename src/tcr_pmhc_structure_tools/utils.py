def get_coords(df):
    return df[['pos_x', 'pos_y', 'pos_z']].to_numpy()
