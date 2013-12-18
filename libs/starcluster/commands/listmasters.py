from completers import ClusterCompleter


class CmdListMaster(ClusterCompleter):
    """
    listmaster [<cluster_tag> ...]

    List master dns address
    """
    names = ['listmaster', 'lm']


    def execute(self, args):
        return self.cm.list_master(cluster_groups=args)
