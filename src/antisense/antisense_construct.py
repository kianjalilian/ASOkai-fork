from abc import ABC, abstractmethod
from sites import TargetSite
from typing import TYPE_CHECKING
from Bio.Seq import Seq
from biochemistry import Chemistry
if TYPE_CHECKING:
    from sites import TargetSite
    
class AntisenseConstruct(ABC):
    """Abstract base class for antisense constructs."""
    def __init__(self, 
                 sequence: Seq,
                 target_site: "TargetSite", 
                 _chemistry: Chemistry,
                 **kwargs):
        """
        Initializes an AntisenseConstruct object.
        
        Args:

        """
        self.id = id

        for key, value in kwargs.items():
            setattr(self, key, value)

    @abstractmethod
    def get_antisense_construct_type(self) -> str:
        pass