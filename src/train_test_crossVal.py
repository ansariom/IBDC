import sys, getopt
import numpy as np
import csv
from sklearn import preprocessing
from numpy import float32, int32
from sklearn.preprocessing.data import StandardScaler
from array import array
import tensorflow as tf

all_features = []
scaled_train_all_features = []
scaled_test_all_features = []
train_all_labels = []
test_all_labels = []
train_all_features = []
train_all_labels = []
batch_index = 0
learning_rate = 0.01
training_epochs = 300
batch_size = 50
display_step = 20
weight_decays = np.arange(0.0001,0.0005,0.0001) 

def next_batch(batch_size, train_features, train_labels):
    global batch_index
    # print ('batch_size', batch_size)
    if batch_index + batch_size > train_labels.shape[0]:
        # print ('resting')
        batch_index = 0
#         train_labels, train_features = shuffle(train_labels, train_features, random_state=np.random.randint(0, 10000))
    batch = (train_labels[batch_index:batch_index + batch_size], train_features[batch_index:batch_index + batch_size])
    batch_index += batch_size
    return batch

def train_test_model(train_features, train_labels, test_features, test_labels, weight_decay):
    feature_size = train_features.shape[1]
    x = tf.placeholder(tf.float32, [None, feature_size])  # mnist data image of shape 28*28=784
    y = tf.placeholder(tf.float32, [None, 1])  # 0-9 digits recognition => 10 classes
    
    'initialize the model weights'
    W = tf.Variable(tf.zeros([feature_size, 1]))
    b = tf.Variable(tf.zeros([1]))
    
    'Construct model'
    logits = (tf.matmul(x, W) + b)  # logit function
    #pred = tf.nn.softmax(logits)   # for more than two class
    pred = tf.nn.sigmoid(logits)
    
    # Minimize error using cross entropy
    #cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y, logits=logits))
    cost = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=y , logits=logits))

    # cost = tf.reduce_mean(-tf.reduce_sum(y * tf.log(pred), reduction_indices=1))
    regularization_loss = tf.nn.l2_loss(W) + tf.nn.l2_loss(b)
    loss = cost + weight_decay * regularization_loss
    # Gradient Descent
    #optimizer = tf.train.GradientDescentOptimizer(learning_rate).minimize(loss)
    optimizer = tf.train.AdamOptimizer().minimize(loss)
    # Initialize the variables (i.e. assign their default value)
    init = tf.global_variables_initializer()
    
    num_examples = train_labels.shape[0]
# Start training
    with tf.Session() as sess:
        coord = tf.train.Coordinator()
        threads = tf.train.start_queue_runners(coord=coord)
        # Run the initializer
        sess.run(init)
        
        # Training cycle
        for epoch in range(training_epochs):
            avg_cost = 0.
            total_batch = int(num_examples / batch_size)
            # Loop over all batches
            for i in range(total_batch):
                (batch_ys, batch_xs) = next_batch(batch_size)
                # Run optimization op (backprop) and cost op (to get loss value)
                _, c = sess.run([optimizer, loss], feed_dict={x: batch_xs,
                                                              y: batch_ys})
                # Compute average loss
                avg_cost += c / total_batch
            # Display logs per epoch step
            if (epoch + 1) % display_step == 0:
                print("Epoch:", '%04d' % (epoch + 1), "cost=", "{:.9f}".format(avg_cost))
    
        print("Optimization Finished!")
    
        
        #print(pred.eval({x: test_features, y: test_labels}))
        # Test model
        #correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
        prediction_op = tf.round(pred)
        correct_prediction = tf.equal(prediction_op, y)
        
        #print(sess.run(correct_prediction, feed_dict={x: test_features, y:test_labels}))
    
        #cm = tf.confusion_matrix(labels=y2, predictions=tf.argmax(pred, 1), num_classes=2)
        # Calculate accuracy
        accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
        print("Param = ", weight_decay, " -- Accuracy:", accuracy.eval({x: test_features, y: test_labels}))


'Split data into train and helout test'
def split_data_train_test(seed):
    print(seed)
    nfeatures = all_features.shape[0]
    ntest = int(0.2 * nfeatures)
    np.random.seed(seed)
    
    test_feature_ids = np.random.choice(nfeatures, ntest, replace = False)
    test_all = all_features[test_feature_ids]
    
    all_idx = np.arange(nfeatures)
    train_ids = np.setdiff1d(all_idx, test_feature_ids)
    train_all = all_features[train_ids]
    
    train_all_labels = np.asanyarray(train_all[:,-1], int32)
    train_all_features = np.asanyarray(train_all[:,1:-1], float32)
    
    test_all_labels = np.asanyarray(test_all[:,-1], int32)
    test_all_features = np.asanyarray(test_all[:,1:-1], float32)
    
    ' scale the train set and use the scaler to scale test set later'
    return train_all_features, train_all_labels

def split_train_cross_val(nfold):
    'split train into a pie with 5 section, 4 train 1 validation'
    ntrain = train_all_features.shape[0]
    pie_size = ntrain/5
    all_idx = np.arange(ntrain)
    
    for i in xrange(nfold):
        start_idx = i * pie_size
        end_idx = start_idx + pie_size
        if (i == nfold - 1):
            end_idx = ntrain
        test_fold_features = train_all_features[start_idx:end_idx]
        test_fold_labels = train_all_labels[start_idx:end_idx]
        
        train_idxs = np.setdiff1d(all_idx, np.arange(start_idx, end_idx))
        train_fold_features = train_all_features[train_idxs]
        train_fold_labels = train_all_labels[train_idxs]
        
        # normalize the data
        scaler = preprocessing.StandardScaler().fit(train_fold_features)
        scaled_train_features = scaler.transform(train_all_features)
        scaled_test_features = scaler.transform(test_fold_features)
        for param in iter(weight_decays):
            print("trying param = ", param)
            train_test_model(scaled_train_features, train_fold_labels, scaled_test_features, test_fold_labels, param)
    
def usage():
    print('train_test_crossval.py -i <input_features.csv> -o <outdir> -n <nfolds> -t <model_type>')
    
def readInputs(argv):
    argsDict = {}
    
    try:
        # Reads Array of tuples such as : [("-i", "in.txt"), ("-o", "out.txt")]
        opts, args = getopt.getopt(argv, "hi:o:n:t:", ["infile=", "outfile=", "nfold=", "modelType="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    if len(argv) < 4:
        usage()
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            usage()  
            sys.exit()
        elif opt in ("-i", "--infile"):
            inputFile = arg
        elif opt in ("-o", "--outfile"):
            outputFile = arg
        elif opt in ("-n", "--nfold"):
            nfold = arg
        elif opt in ("-t", "--modelType"):
            model_type = arg
            
    argsDict["infile"] = inputFile
    argsDict["outfile"] = outputFile
    argsDict["nfold"] = nfold
    argsDict["model_type"] = model_type
    return argsDict

if __name__ == "__main__":
    first_line = True
    argsDict = readInputs(sys.argv[1:])
    with open(argsDict["infile"], 'rt') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            if (first_line):
                first_line = False
                continue
            all_features.append(np.asarray(row))
        all_features = np.asanyarray(all_features)

    seed = np.random.randint(12000)
    train_all_features, train_all_labels = split_data_train_test(seed)
    split_train_cross_val(5)
    