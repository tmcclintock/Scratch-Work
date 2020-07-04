import torch

device = torch.device('cpu') #'cuda' if we have a GPU

N, D_in, H, D_out = 64, 1000, 100, 10

x = torch.randn(N, D_in, device = device)
y = torch.randn(N, D_out, device = device)

model = torch.nn.Sequential(
    torch.nn.Linear(D_in, H),
    torch.nn.ReLU(),
    torch.nn.Linear(H, D_out),
).to(device)

loss_fn = torch.nn.MSELoss(reduction = 'sum')

learning_rate = 1e-4
optimizer = torch.optim.Adam(model.parameters(), lr = learning_rate)
for t in range(500):
    y_pred = model(x)

    loss = loss_fn(y_pred, y)
    print(f"Epoch {t} loss = {loss.item():.3e}")

    """
    model.zero_grad()
    loss.backward()
    with torch.no_grad():
        for param in model.parameters():
            param.data -= learning_rate * param.grad
    """

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()