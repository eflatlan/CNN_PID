import tensorflow as tf
from tensorflow.keras.layers import Input, Conv2D, Flatten, Dense, concatenate
from tensorflow.keras.models import Model

# Pion, Kaon, Proton, Unkown
num_classes = 4

# Define input shapes
map_shape = (144, 144, 1) 
momentum_shape = (1,)
azimuthal_angle_shape = (1,)
polar_angle_shape = (1,)

# Define inputs
map_input = Input(shape=map_shape, name='map_input')
momentum_input = Input(shape=momentum_shape, name='momentum_input')
azimuthal_angle_input= Input(shape=azimuthal_angle_shape, name='azimuthal_angle_shape')
polar_angle_input = Input(shape=polar_angle_shape, name='polar_angle_shape')

# Define convolutional layers for each map input
conv1 = Conv2D(32, (3, 3), activation='relu')(map_input)
conv2 = Conv2D(32, (3, 3), activation='relu')(conv1)
conv3 = Conv2D(32, (3, 3), activation='relu')(conv2)
flat_map = Flatten()(conv3)

# Concatenate map features with other inputs
concat = concatenate([flat_map, momentum_input, azimuthal_angle_input, polar_angle_input])

# Define fully connected layers
fc1 = Dense(128, activation='relu')(concat)
fc2 = Dense(64, activation='relu')(fc1)
output = Dense(num_classes, activation='softmax')(fc2) 

# Define the model
model = Model(inputs=[map_input, momentum_input, azimuthal_angle_input, azimuthal_angle_input], outputs=output)

# Compile the model
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
