import tensorflow as tf
from tensorflow.keras.layers import Input, Conv2D, Flatten, Dense, concatenate, Lambda
from tensorflow.keras.models import Model
import tensorflow.keras.backend as K

# Pion, Kaon, Proton, Unkown
num_classes = 4

# Define input shapes
map_shape = (160, 144, 7)  # 160 x pads [0..159] | 144 y pads [0..143] | 7 chambers 
momentum_shape = (1,) # p_T
azimuthal_angle_shape = (1,) # theta
polar_angle_shape = (1,) # phi_p
MIP_impact_shape = (2,) # impact point of MIP {xP, yP}
ref_index_shape = (1,) # refractive index (n)
# + Track inclination?

# Define inputs
map_input = Input(shape=map_shape, name='map_input')
momentum_input = Input(shape=momentum_shape, name='momentum_input')
azimuthal_angle_input= Input(shape=azimuthal_angle_shape, name='azimuthal_angle_input')
polar_angle_input = Input(shape=polar_angle_shape, name='polar_angle_input')
MIP_impact_input = Input(shape=MIP_impact_shape, name='MIP_impact_input')
ref_index_input = Input(shape=ref_index_shape, name='refractive_index_input')

'''# Define input layers for position
pos_x_input = Input(shape=(1,), name='pos_x_input')
pos_y_input = Input(shape=(1,), name='pos_y_input')

# Define a lambda layer to convert position inputs to one-hot encoded vectors
def one_hot(x):
    return K.one_hot(K.cast(x, dtype='int32'), num_classes=144)

# Apply the one-hot encoding to position inputs
# Apply the one-hot encoding to position inputs
one_hot_x = Lambda(lambda x: K.one_hot(K.cast(x, dtype='int32'), num_classes))(pos_x_input)
one_hot_y = Lambda(lambda x: K.one_hot(K.cast(x, dtype='int32'), num_classes))(pos_y_input)
one_hot_x = Lambda(lambda x: K.expand_dims(x, axis=1))(one_hot_x)
one_hot_y = Lambda(lambda x: K.expand_dims(x, axis=1))(one_hot_y)'''

track_inputs = [momentum_input, azimuthal_angle_input, polar_angle_input, ref_index_input] 
#track_inputs = [momentum_input, azimuthal_angle_input, polar_angle_input, one_hot_x, one_hot_y, ref_index_input] 

# Define convolutional layers for each map input
conv1 = Conv2D(32, (3, 3), activation='relu')(map_input)
conv2 = Conv2D(32, (3, 3), activation='relu')(conv1)
conv3 = Conv2D(32, (3, 3), activation='relu')(conv2)
flat_map = Flatten()(conv3)

# Concatenate map features with other inputs
#concat = concatenate([flat_map, momentum_input, azimuthal_angle_input, polar_angle_input, one_hot_x, one_hot_y, ref_index_input])
concat = concatenate([flat_map, momentum_input, azimuthal_angle_input, polar_angle_input, ref_index_input, MIP_impact_input])

# Define fully connected layers
fc1 = Dense(128, activation='relu')(concat)
fc2 = Dense(64, activation='relu')(fc1)
output = Dense(num_classes, activation='softmax')(fc2) 

# Define the model
model = Model(inputs=[map_input, momentum_input, azimuthal_angle_input, polar_angle_input, ref_index_input, MIP_impact_input], outputs=output)

# Compile the model
model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
