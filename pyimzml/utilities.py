"""Utility functions"""


def chunks(item_list, n_items):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(item_list), n_items):
        yield item_list[i : i + n_items]
