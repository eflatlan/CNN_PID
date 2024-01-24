from tensorflow.keras.callbacks import Callback

class GradualBatchSizeIncrease(Callback):
    def __init__(self,  start_size=16, max_size=256, increase_factor=2, interval=20):
        super(GradualBatchSizeIncrease, self).__init__()
        self.batch_size = start_size
        self.max_size = max_size
        self.increase_factor = increase_factor
        self.interval = interval

    def on_epoch_end(self,  epoch, logs=None):
        if (epoch + 1) % self.interval == 0 and self.batch_size < self.max_size:
            self.batch_size = min(self.batch_size * self.increase_factor, self.max_size)
            self.model.batch_size = self.batch_size
            print(f"\nEpoch {epoch + 1}: Increasing batch size to {self.batch_size}.\n")
