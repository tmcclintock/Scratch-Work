#A two layer NN in just numpy
import numpy as np

N, D_in, H, D_out = 64, 1000, 100, 10

x = np.random.randn(N, D_in)
y = np.random.randn(N, D_out)

w1 = np.random.randn(D_in, H)
w2 = np.random.randn(H, D_out)

learning_rate = 1e-6
for t in range(500):
    #Forward pass: compute predicted targets
    h = x.dot(w1) #(N, H)
    h_relu = np.maximum(h, 0) #(N, H)
    y_pred = h_relu.dot(w2) #(N, D_out)

    #Compute and print loss
    loss = np.square(y_pred - y).sum() # scalar
    #print(f"Epoch {t} loss = {loss:.3e}")

    #Back prop
    grad_y_pred = 2 * (y_pred - y) #(N, D_out)
    grad_w2 = h_relu.T.dot(grad_y_pred) #(H, D_out)
    grad_h_relu = grad_y_pred.dot(w2.T) #(N, H)
    grad_h = grad_h_relu.copy()
    grad_h[h < 0] = 0 #(N, H)
    grad_w1 = x.T.dot(grad_h) #(D_in, H)

    #Update
    w1 -= learning_rate * grad_w1
    w2 -= learning_rate * grad_w2
