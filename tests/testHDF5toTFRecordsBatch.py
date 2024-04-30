import tensorflow as tf
import h5py
import numpy as np

def _bytes_feature(value):
    """Returns a bytes_list from a string / byte."""
    if isinstance(value, type(tf.constant(0))):
        value = value.numpy() # BytesList won't unpack a string from an EagerTensor.
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def serialize_example(real_field, imag_field, tbp):
    """
    Creates a tf.Example message ready to be written to a file.
    """
    # Create a dictionary mapping the feature name to the tf.Example-compatible
    # data type.
    feature = {
        'tbp': _bytes_feature(tf.io.serialize_tensor(tbp)),
        'real_field': _bytes_feature(tf.io.serialize_tensor(real_field)),
        'imag_field': _bytes_feature(tf.io.serialize_tensor(imag_field)),
    }

    # Create a Features message using tf.train.Example.
    example_proto = tf.train.Example(features=tf.train.Features(feature=feature))
    return example_proto.SerializeToString()

if __name__ == '__main__':
    N = 128
    NUMBER_OF_PULSES = 20000
    BATCH_SIZE = 1024

    nBatches = np.ceil(NUMBER_OF_PULSES / BATCH_SIZE).astype(int)

    for batch in range(nBatches):
        filename = f'{NUMBER_OF_PULSES}_randomNormalizedPulses_N{N}_batch{batch}'
        with h5py.File(filename + '.h5', 'r') as f:
            with tf.io.TFRecordWriter(filename + '.tfrecords') as writer:
                for group_name in f:
                    group = f[group_name]
                    
                    # Read TBP
                    TBP = group.attrs['TBP']
                    
                    # Read real part of field
                    real_field = np.array(group['real_field'])
                    
                    # Read imaginary part of field
                    imag_field = np.array(group['imag_field'])
                    
                    example = serialize_example(real_field, imag_field, TBP)
                    writer.write(example)