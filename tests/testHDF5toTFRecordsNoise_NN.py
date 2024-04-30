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


def serialize_example(real_field, imag_field, noisy_trace, real_retrieved_field, imag_retrieved_field):
    """
    Creates a tf.Example message ready to be written to a file.
    """
    # Create a dictionary mapping the feature name to the tf.Example-compatible
    # data type.
    feature = {
        'real_field': _bytes_feature(tf.io.serialize_tensor(real_field)),
        'imag_field': _bytes_feature(tf.io.serialize_tensor(imag_field)),
        'noisy_trace': _bytes_feature(tf.io.serialize_tensor(noisy_trace)),
        'real_retrieved_field': _bytes_feature(tf.io.serialize_tensor(real_retrieved_field)),
        'imag_retrieved_field': _bytes_feature(tf.io.serialize_tensor(imag_retrieved_field))

    }

    # Create a Features message using tf.train.Example.
    example_proto = tf.train.Example(
        features=tf.train.Features(feature=feature))
    return example_proto.SerializeToString()

def read_tfrecord_noisyTraces(FILE_PATH, N, NUMBER_OF_PULSES, BATCH_SIZE):
    """
    Read the TFRecord file and return a dataset with the pulses and their SHG-FROG noisy traces.
    The pulses come already normalized.
    Traces are computed and normalized in this function to compare to the noisy traces.


    Args:
        FILE_PATH: str
            Path to the TFRecord file
        N: int
            Number of points in the SHG-FROG trace
        NUMBER_OF_PULSES: int
            Number of pulses in the database
        BATCH_SIZE: int
            Size of the batches to use in the dataset

    Returns:
        dataset: tf.data.Dataset
            Dataset with the pulse database
    """
    # Create a dataset from the TFRecord file
    raw_dataset = tf.data.TFRecordDataset(FILE_PATH)

    # Map the parse function over the dataset
    parsed_dataset = raw_dataset.map(_parse_function)

    # Create empty datasets
    field_dataset = tf.data.Dataset.from_tensor_slices([])
    noisy_trace_dataset = tf.data.Dataset.from_tensor_slices([])
    retrieved_field_dataset = tf.data.Dataset.from_tensor_slices([])

    for real_field, imag_field, noisy_trace in parsed_dataset:
        # Concat the real and imaginary parts into a single tensor
        pulse = tf.concat([tf.reshape(real_field, (1, N)), tf.reshape(imag_field, (1,N))], axis=1)
        
        # Add the pulse to the field_dataset
        field_dataset = field_dataset.concatenate(
            tf.data.Dataset.from_tensor_slices(tf.cast(pulse, tf.float32)))

        # Concat the real and imaginary parts into a single tensor
        retrieved_field = tf.concat([tf.reshape(real_retrieved_field, (1, N)), tf.reshape(imag_retrieved_field, (1,N))], axis=1)

        # Add the retrieved field to the retrieved_field_dataset
        retrieved_field_dataset = retrieved_field_dataset.concatenate(
            tf.data.Dataset.from_tensor_slices(tf.cast(retrieved_field, tf.float32)))

        # Add the noisy trace to the noisy_trace_dataset
        noisy_trace_dataset = noisy_trace_dataset.concatenate(
            tf.data.Dataset.from_tensor_slices(tf.cast(noisy_trace, tf.float32)))


    # Compute trace of the electric fields
    fourier = fourier_utils(N, 1/N)
    # We pass it as a batch for performance (it's faster to compute the fourier transform of a batch of pulses)
    # and less memory intensive
    trace_dataset = field_dataset.batch(BATCH_SIZE).map(fourier.compute_trace)

    # Norm the traces by dividing by the maximum value of each trace
    trace_dataset = trace_dataset.map(lambda x: x / tf.reduce_max(tf.abs(x), axis=[1, 2], keepdims=True))

    # Unbatch the dataset
    trace_dataset = trace_dataset.unbatch()

    dataset = tf.data.Dataset.zip((noisy_trace_dataset, trace_dataset, field_dataset, retrieved_field_dataset))

    return dataset


if __name__ == '__main__':
    N = 128
    NUMBER_OF_PULSES = 10
    SNR = 10

    filename = f'{NUMBER_OF_PULSES}_randomPulses_N{N}_{SNR}SNR'
    with h5py.File(filename + '.h5', 'r') as f:
        with tf.io.TFRecordWriter(filename + '.tfrecords') as writer:
            for key in f.keys():

                real_original_field = f[key]['real_original_field'][:]

                imag_original_field = f[key]['imag_original_field'][:]

                original_trace = f[key]['original_trace'][:]

                noisy_trace = f[key]['noisy_trace'][:]

                # # Subplot with original and noisy traces imshow
                # fig, axs = plt.subplots(1, 2)
                # # fig.suptitle('Original and Noisy Traces')
                # axs[0].imshow(original_trace, cmap='nipy_spectral')
                # # Unset xticks and yticks
                # axs[0].set_yticks([])
                # axs[0].set_xticks([])
                # axs[0].set_xlabel('Delays')
                # axs[0].set_ylabel('Frequencies')
                # axs[1].imshow(noisy_trace, cmap='nipy_spectral')
                # # Unset xticks and yticks
                # axs[1].set_yticks([])
                # axs[1].set_xticks([])
                # axs[1].set_xlabel('Delays')
                # axs[1].set_ylabel('Frequencies')

                # plt.show()

                real_retrieved_field = f[key]['real_retrieved_field'][:]

                imag_retrieved_field = f[key]['imag_retrieved_field'][:]

                example = serialize_example(real_original_field, imag_original_field, noisy_trace, real_retrieved_field, imag_retrieved_field)

                # Write the `tf.Example` observations to the file.
                writer.write(example)

