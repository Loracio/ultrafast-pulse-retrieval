import tensorflow as tf
import h5py
import numpy as np
import matplotlib.pyplot as plt
from plotRetrievalResults import *
import scienceplots
plt.style.use(['science', 'ieee'])


def _bytes_feature(value):
    """Returns a bytes_list from a string / byte."""
    if isinstance(value, type(tf.constant(0))):
        # BytesList won't unpack a string from an EagerTensor.
        value = value.numpy()
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))


def serialize_example(real_field, imag_field, original_trace, real_retrieved_field, imag_retrieved_field, noisy_trace):
    """
    Creates a tf.Example message ready to be written to a file.
    """
    # Create a dictionary mapping the feature name to the tf.Example-compatible
    # data type.
    feature = {
        'real_field': _bytes_feature(tf.io.serialize_tensor(real_field)),
        'imag_field': _bytes_feature(tf.io.serialize_tensor(imag_field)),
        'original_trace': _bytes_feature(tf.io.serialize_tensor(original_trace)),
        'real_retrieved_field': _bytes_feature(tf.io.serialize_tensor(real_retrieved_field)),
        'imag_retrieved_field': _bytes_feature(tf.io.serialize_tensor(imag_retrieved_field)),
        'noisy_trace': _bytes_feature(tf.io.serialize_tensor(noisy_trace)),
    }

    # Create a Features message using tf.train.Example.
    example_proto = tf.train.Example(
        features=tf.train.Features(feature=feature))
    return example_proto.SerializeToString()


if __name__ == '__main__':
    N = 128
    NUMBER_OF_PULSES = 10

    filename = f'{NUMBER_OF_PULSES}_randomPulses_N{N}with_0.000000noise'
    with h5py.File(filename + '.h5', 'r') as f:
        with tf.io.TFRecordWriter(filename + '.tfrecords') as writer:
            for key in f.keys():

                real_original_field = f[key]['real_original_field'][:]

                imag_original_field = f[key]['imag_original_field'][:]

                original_trace = f[key]['original_trace'][:]

                noisy_trace = f[key]['noisy_trace'][:]

                # Subplot with original and noisy traces imshow
                fig, axs = plt.subplots(1, 2)
                # fig.suptitle('Original and Noisy Traces')
                axs[0].imshow(original_trace, cmap='nipy_spectral')
                # Unset xticks and yticks
                axs[0].set_yticks([])
                axs[0].set_xticks([])
                axs[0].set_xlabel('Delays')
                axs[0].set_ylabel('Frequencies')
                axs[1].imshow(noisy_trace, cmap='nipy_spectral')
                # Unset xticks and yticks
                axs[1].set_yticks([])
                axs[1].set_xticks([])
                axs[1].set_xlabel('Delays')
                axs[1].set_ylabel('Frequencies')

                plt.show()

                real_retrieved_field = f[key]['real_retrieved_field'][:]

                imag_retrieved_field = f[key]['imag_retrieved_field'][:]

                example = serialize_example(real_original_field, imag_original_field, original_trace, real_retrieved_field, imag_retrieved_field, noisy_trace)

                # Write the `tf.Example` observations to the file.
                writer.write(example)

