seq = ("TAT"
       "ACGACGACG"
       "GGT"
       "ACGACGACG"
       "GGTGGTGGT"
       "GTAGTGGGTGAGTGCGTGC")


class CustomQueue:
    def __init__(self, k):
        self.k = k
        queue_size = 1
        self.queues = [["$"] * queue_size for _ in range(k)]
        self.frame = 0
        self.scores = [[] * queue_size for _ in range(k)]

    def enqueue(self, kmer):

        queue = self.queues[self.frame]
        # print("frame", self.frame, "kmer", kmer, "queue", queue, "score", kmer in queue)
        if kmer in queue:
            self.scores[self.frame].append(1)
        else:
            self.scores[self.frame].append(0)

        if len(queue) >= 1:  # and kmer not in queue:
            queue.pop(0)
            queue.append(kmer)
        self.frame = (self.frame + 1) % self.k

    def get_repeat(self, sequence):
        for frame in range(self.k):
            result = ""
            last_kmer = ""
            for i in range(len(self.scores[frame])):
                position_in_seq = i * self.k + frame
                current_kmer = sequence[position_in_seq: position_in_seq + self.k]
                next_kmer = sequence[position_in_seq + self.k: position_in_seq + self.k * 2]

                if self.scores[frame][i] == 0:
                    if current_kmer == last_kmer:
                        last_kmer = ""
                    else:
                        result += current_kmer[::-1]
                if self.scores[frame][i] == 1:
                    if last_kmer == current_kmer:
                        result = result[:-2 - k - 2] + " " + str(int(result[-4 - k]) + 1)
                    else:
                        result += " 2_)" + current_kmer[::-1] + "("

                last_kmer = current_kmer

            print("k =", k, "frame =", frame)
            print(result[::-1])


kmers_to_check = [2, 3, 4, 5]
custom_queues = {k: CustomQueue(k) for k in kmers_to_check}
max_k = max(kmers_to_check)
last_max_kmer = "$" * max_k
for base in seq:
    last_max_kmer = last_max_kmer[1:] + base
    for k in kmers_to_check:
        current_kmer = last_max_kmer[-k:]
        if "$" not in current_kmer:
            custom_queues[k].enqueue(current_kmer)

for k, queue in custom_queues.items():
    print(seq)
    # print(" " * (k - 2), "".join([str(i) for i in queue.score]))
    # for i in range(len(queue.score) - k - k + 1):
    #     if queue.score[i + k + k - 1] != 0 and queue.score[i] == 0:
    #         queue.score[i:i + k + k - 1] = [queue.score[i + k + k - 1]] * (k + k - 1)
    queue.get_repeat(seq)
    # print("####")
    # mask = "".join([" " if i == 0 else str(k) for i in queue.score])
    # print(" " * (k - 2), mask)
    # print()
