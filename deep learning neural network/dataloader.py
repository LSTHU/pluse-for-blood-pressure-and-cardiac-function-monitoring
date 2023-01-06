import torch
import pandas as pd
import numpy as np
import random
from torch.utils.data import Dataset, DataLoader, random_split, TensorDataset
from sklearn.preprocessing import StandardScaler

length = 32
seed = 12345

torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
random.seed(seed)


class Custom_Dataset(Dataset):
    def __init__(self, mode='train'):
        super(Custom_Dataset, self).__init__() 
        data = pd.read_csv('pulse.csv', header=None)
        data[data.columns[9:]] = StandardScaler().fit_transform(data[data.columns[9:]])
        self.data = data.values
    def __len__(self):
        return len(self.data) - length

    def __getitem__(self, item):
        idx = item % (len(self.data) - length)
        x = self.data[idx: idx+length, 9:]
        y = self.data[idx+length, 7:9]

        return torch.unsqueeze(torch.tensor(x).float(), 0), torch.tensor(y).float()

def get_dataloader():
    dataset = Custom_Dataset()
#
    train_dataset, valid_dataset = random_split(dataset, (int(0.8 * len(dataset)), len(dataset) - int(0.8 * len(dataset))))
    train_dl = DataLoader(train_dataset, batch_size=16, num_workers=4, shuffle=True)
    valid_dl = DataLoader(valid_dataset, batch_size=16, num_workers=4, shuffle=False)
    return train_dl, valid_dl

if __name__ == '__main__':
    train_dl, valid_dl = get_dataloader()
    for item in valid_dl:
        print(item[1].shape)
