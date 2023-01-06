import time
import torch
import numpy as np
import pandas as pd
import torch.nn as nn
import random
from sklearn.metrics import mean_absolute_error, r2_score
from dataloader0720 import get_dataloader, Custom_Dataset
import matplotlib.pyplot as plt

epochs = 200
seed = 12345

torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
random.seed(seed)


class Block(nn.Module):
    def __init__(self, in_channel, out_channel):
        super(Block, self).__init__()
        self.module = nn.Sequential(
            nn.Conv2d(in_channel, out_channel, kernel_size=(3, 3), padding=(1, 1)),
            nn.LeakyReLU(),
            nn.BatchNorm2d(out_channel),
        )

    def forward(self, x):
        return self.module(x)


class Network(nn.Module):
    def __init__(self):
        super(Network, self).__init__()
        self.block1 = Block(1, 64)
        self.block2 = Block(64, 128)
        self.block3 = Block(128, 256)
        self.block4 = Block(256, 256)

        self.dropout = nn.Dropout(0.3)
        self.maxpool = nn.MaxPool2d((2, 8), (2, 8))
        self.clf = nn.Sequential(
            nn.Linear(2048, 512),
            nn.ReLU(),
            nn.BatchNorm1d(512),
            nn.Dropout(0.4),
            nn.Linear(512, 2),
        )

    def forward(self, x):
        x = self.block1(x)
        x = self.block2(x)
        x = self.maxpool(x)
        x = self.block3(x)
        x = self.maxpool(x)
        x = self.block4(x)
        x = self.maxpool(x)

        x = torch.flatten(x, 1)
        x = self.dropout(x)
        x = self.clf(x)
        return x


def train(train_dl, valid_dl):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = Network().to(device)
    loss_func = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    best_loss = -1000
    history_loss = [] # history
    for epoch in range(epochs):
        start = time.time()
        model.train()
        preds, labels = [], []
        for item in train_dl:
            optimizer.zero_grad()
            data, label = item[0].to(device), item[1].to(device)
            pred = model(data) # output
            loss = loss_func(pred, label) # loss
            loss.backward() # BP
            optimizer.step() # update weight
            preds.extend(model(data).detach().cpu().numpy())
            labels.extend(label.detach().cpu().numpy())

        model.eval()
        test_preds, test_labels = [], []
        for item in valid_dl:
            data, label = item[0].to(device), item[1].to(device)
            test_preds.extend(model(data).detach().cpu().numpy())
            test_labels.extend(label.detach().cpu().numpy())

        print('MAE: {:.4f}, R2: {:.4f}, MAE: {:.4f}, R2: {:.4f}'.format(mean_absolute_error(labels, preds),
                                                                        r2_score(labels, preds),
                                                                        mean_absolute_error(test_labels, test_preds),
                                                                        r2_score(test_labels, test_preds)))
        history_loss.append(mean_absolute_error(labels, preds))
        if best_loss < r2_score(test_labels, test_preds):
            print(r2_score(test_labels, test_preds))
            best_loss = r2_score(test_labels, test_preds)
            torch.save(model, 'BP_model.pth')

        end = time.time()
        print(f"one epoch time {end - start} s")

    model = torch.load('BP_model.pth',
                       map_location=torch.device(device))
    model.eval()
    test_preds, test_labels = [], []
    losses = []
    for item in valid_dl:
        data, label = item[0].to(device), item[1].to(device)
        test_preds.extend(model(data).detach().cpu().numpy())
        test_labels.extend(label.detach().cpu().numpy())
    print()
    print('MAE: {:.4f}, R2: {:.4f}, MAE: {:.4f}, R2: {:.4f}'.format(mean_absolute_error(labels, preds),
                                                                    r2_score(labels, preds),
                                                                    mean_absolute_error(test_labels, test_preds),
                                                                    r2_score(test_labels, test_preds)))

    plt.figure()
    plt.plot(np.arange(epochs) + 1, history_loss, label='loss')

    pd.DataFrame(columns=['loss'], data=history_loss).to_csv('BP_loss.csv', index=False)
    plt.savefig('BP_loss.png')


def predict(valid_dl):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = torch.load('BP_model.pth',
                       map_location=torch.device(device))
    model.eval()
    test_preds, test_labels = [], []
    losses = []
    for item in valid_dl:
        data, label = item[0].to(device), item[1].to(device)
        test_preds.extend(model(data).detach().cpu().numpy())
        test_labels.extend(label.detach().cpu().numpy())

    print('MAE: {:.4f}, R2: {:.4f}'.format(mean_absolute_error(test_labels, test_preds),
                                           r2_score(test_labels, test_preds)))
    test_preds, test_labels = [], []
    dataset = Custom_Dataset()
    dl = torch.utils.data.DataLoader(dataset, batch_size=32, shuffle=False)
    for item in dl:
        data, label = item[0].to(device), item[1].to(device)
        test_preds.extend(model(data).detach().cpu().numpy())
        test_labels.extend(label.detach().cpu().numpy())

    pd.DataFrame(columns=['prediction1', 'prediction2'], data=test_preds).to_csv('BP_prediction.csv', index=False)
    pd.DataFrame(columns=['label1', 'label2'], data=test_labels).to_csv('BP_label.csv', index=False)


if __name__ == '__main__':
    train_dl, valid_dl = get_dataloader()
    #train(train_dl, valid_dl)
    predict(valid_dl)
