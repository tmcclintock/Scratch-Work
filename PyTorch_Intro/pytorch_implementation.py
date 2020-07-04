import torch

class MyReLU(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x):
        ctx.save_for_backward(x)
        return x.clamp(min = 0)

    def backward(ctx, grad_output):
        x, = ctx.saved_tensors
        grad_x = grad_output.clone()
        grad_x[x < 0] = 0
        return grad_x

device = torch.device('cpu')
#device = torch.device('cuda') #if we had a GPU

N, D_in, H, D_out = 64, 1000, 100, 10

x = torch.randn(N, D_in, device = device)
y = torch.randn(N, D_out, device = device)

w1 = torch.randn(D_in, H, device = device, requires_grad = True)
w2 = torch.randn(H, D_out, device = device, requires_grad = True)

learning_rate = 1e-6
for t in range(500):
    #Forward pass
    #h = x.mm(w1) #(N, H)
    #h_relu = h.clamp(min = 0) #(N, H)
    #y_pred = h_relu.mm(w2) #(N, D_out)
    y_pred = MyReLU.apply(x.mm(w1)).mm(w2)
    
    #Compute and print loss
    #note that loss is a scalar with shape ()
    #we get it's value as a float with loss.item()
    loss = (y_pred - y).pow(2).sum()
    print(f"Epoch {t} loss = {loss.item():.3e}")

    #Backprop manually for now
    """
    grad_y_pred = 2.0 * (y_pred - y) #(N, D_out)
    grad_w2 = h_relu.t().mm(grad_y_pred) #(H, D_out)
    grad_h_relu = grad_y_pred.mm(w2.t()) #(N, H)
    grad_h = grad_h_relu.clone()
    grad_h[h < 0] = 0 #(N, H)
    grad_w1 = x.t().mm(grad_h) #(D_in, H)

    #Update weights
    w1 -= learning_rate * grad_w1
    w2 -= learning_rate * grad_w2
    """
    #Backprop w/ autograd
    loss.backward()
    with torch.no_grad():
        w1 -= learning_rate * w1.grad
        w2 -= learning_rate * w2.grad

        w1.grad.zero_()
        w2.grad.zero_()
